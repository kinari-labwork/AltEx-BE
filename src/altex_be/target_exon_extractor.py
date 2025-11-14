from __future__ import annotations
import pandas as pd
import uuid
import logging
from . import logging_config # noqa: F401

# BED形式も0base-start, 1base-endであるため、refFlatのexonStartsとexonEndsをそのまま使用する

def explode_classified_refflat(classified_refflat: pd.DataFrame, target_exon: str = "all") -> pd.DataFrame:
    classified_refflat = classified_refflat.explode(["exonStarts", "exonEnds", "exontype", "exon_position", "exonlengths", "flame"])
    classified_refflat["cds_info"] = classified_refflat["cds_info"].apply(lambda x: set(x))  # set型に変換
    classified_refflat = classified_refflat.drop(columns = ["exons"])
    classified_refflat[["exonStarts", "exonEnds"]] = classified_refflat[["exonStarts", "exonEnds"]].astype(
        int
    )  # int型に変換
    # exontypeがskippedまたはuniqueのエキソンだけを抽出
    if target_exon == "alternative_exons":
        classified_refflat = classified_refflat[classified_refflat["exontype"].apply(lambda x: x in ("skipped", "unique","a3ss-long","a5ss-long"))]
    elif target_exon == "all":
        # exon数が2以下かつ、すべてのエキソンがconstitutiveである遺伝子は、スプライシング操作をしても意味がないため除外する
        def filter_genes(group):
            num_exons = len(group)
            has_alternative = any(group['exontype'] != 'constitutive')
            if num_exons <= 2 and not has_alternative:
                return False
            return True
        classified_refflat = classified_refflat.groupby('geneName').filter(filter_genes)
        classified_refflat = classified_refflat[classified_refflat["exontype"].apply(lambda x: x in ("skipped", "unique","a3ss-long","a5ss-long","constitutive"))]
    # 重複を削除し一方だけ残す
    classified_refflat = classified_refflat.drop_duplicates(subset=["chrom", "exonStarts", "exonEnds"])
    classified_refflat['uuid'] = [uuid.uuid4().hex for _ in range(len(classified_refflat))]  # 一意のIDを生成
    return classified_refflat.reset_index(drop=True)

def format_classified_refflat_to_bed(exploded_classified_refflat: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        スプライシングイベントに応じてアノテーションしたrefFlatのデータフレームから
        編集対象となるエキソンだけを抽出し、1エキソン1行のデータフレームに変換する
    Parameters:
        data: pd.DataFrame, exontypeをアノテーション済みのrefFlatのデータフレーム
    Returns:
        pd.DataFrame
    """
    exploded_classified_refflat = exploded_classified_refflat[
        [
            "chrom",
            "strand",
            "exonStarts",
            "exonEnds",
            "exontype",
            "exon_position",
            "uuid"
        ]
    ].copy()
    # 編集のために、リストになっている列を展開する
    exploded_classified_refflat.loc[:, 'score'] = 0  # BED形式のスコア列を追加
    #BED に合わせたカラム順に並べ替え
    exploded_classified_refflat = exploded_classified_refflat.rename(columns={"uuid": "name"})
    classified_refflat = exploded_classified_refflat[["chrom", "exonStarts", "exonEnds", "name", "score", "strand", "exontype", "exon_position"]]
    return classified_refflat.reset_index(drop=True)

def extract_splice_acceptor_regions(target_exon_df: pd.DataFrame, window: int) -> pd.DataFrame:
    """
    Purpose :
        抜き出したskipped or unique exonのexonStart/Endから、SA部位周辺の、windowで指定した幅の座位を示すDataFrameを作成する
        strandが+の時はexonStartがSplice Acceptor, -の時はその逆でexonEndがSAになる
    """
    splice_acceptor_single_exon_df = target_exon_df.copy()
    splice_acceptor_single_exon_df["chromStart"] = splice_acceptor_single_exon_df.apply(
        lambda row: row["exonStarts"] - window
        if row["strand"] == "+"
        else row["exonEnds"] - window,
        axis=1,
    )
    splice_acceptor_single_exon_df["chromEnd"] = splice_acceptor_single_exon_df.apply(
        lambda row: row["exonStarts"] + window
        if row["strand"] == "+"
        else row["exonEnds"] + window,
        axis=1,
    )
    return splice_acceptor_single_exon_df[["chrom","chromStart","chromEnd","name","score","strand"]].reset_index(drop=True)


def extract_splice_donor_regions(target_exon_df: pd.DataFrame, window: int) -> pd.DataFrame:
    """
    Purpose :
        抜き出したskipped or unique exonのexonStart/Endから、SD部位周辺の、windowで指定した幅の座位を示すDataFrameを作成する
        strandが+の時はexonEndがSplice Donor, -の時はその逆でexonStartがSDになる
    """
    splice_donor_single_exon_df = target_exon_df.copy()
    splice_donor_single_exon_df["chromStart"] = splice_donor_single_exon_df.apply(
        lambda row: row["exonEnds"] - window
        if row["strand"] == "+"
        else row["exonStarts"] - window,
        axis=1,
    )
    splice_donor_single_exon_df["chromEnd"] = splice_donor_single_exon_df.apply(
        lambda row: row["exonEnds"] + window
        if row["strand"] == "+"
        else row["exonStarts"] + window,
        axis=1,
    )
    return splice_donor_single_exon_df[["chrom","chromStart","chromEnd","name","score","strand"]].reset_index(drop=True)

def wrap_extract_target_exon(classified_refflat: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Purpose:
    このモジュールの操作をまとめて実行するためのラッパー関数
    """
    exploded_classified_refflat = explode_classified_refflat(classified_refflat, target_exon="all")
    target_exon_df = format_classified_refflat_to_bed(exploded_classified_refflat)
    if target_exon_df is None:
        logging.warning("there are no exons in your interested genes which have at least one targetable splicing event")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    splice_acceptor_single_exon_df = extract_splice_acceptor_regions(target_exon_df, 25)
    splice_donor_single_exon_df = extract_splice_donor_regions(target_exon_df, 25)
    return splice_acceptor_single_exon_df, splice_donor_single_exon_df, exploded_classified_refflat