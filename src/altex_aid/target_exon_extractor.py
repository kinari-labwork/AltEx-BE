from __future__ import annotations

import pandas as pd

import uuid

# BED形式も0base-start, 1base-endであるため、refFlatのexonStartsとexonEndsをそのまま使用する


def extract_target_exon(classified_refflat: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        スプライシングイベントに応じてアノテーションしたrefFlatのデータフレームから
        編集対象となるエキソンだけを抽出し、1エキソン1行のデータフレームに変換する
    Parameters:
        data: pd.DataFrame, exontypeをアノテーション済みのrefFlatのデータフレーム
    Returns:
        pd.DataFrame
    """
    # explodeで増加する行数を抑えるために、先に少なくともskipped exonまたはunique exonを1つ以上持つトランスクリプトを抽出
    classified_refflat = classified_refflat[
        [
            "chrom",
            "strand",
            "exonStarts",
            "exonEnds",
            "exontype",
            "exon_position",
        ]
    ]
    # 編集のために、リストになっている列を展開する
    classified_refflat = classified_refflat.explode(["exonStarts", "exonEnds", "exontype", "exon_position"])
    classified_refflat[["exonStarts", "exonEnds"]] = classified_refflat[["exonStarts", "exonEnds"]].astype(
        int
    )  # int型に変換
    # exontypeがskippedまたはuniqueのエキソンだけを抽出
    classified_refflat = classified_refflat[classified_refflat["exontype"].apply(lambda x: x in ("skipped", "unique","a3ss-long","a5ss-long"))]
    # 重複を削除し一方だけ残す
    classified_refflat = classified_refflat.drop_duplicates(subset=["chrom", "exonStarts", "exonEnds"])
    classified_refflat['name'] = [uuid.uuid4().hex for _ in range(len(classified_refflat))]  # 一意のIDを生成
    classified_refflat['score'] = 0  # BED形式のスコア列を追加
    #BED に合わせたカラム順に並べ替え
    classified_refflat = classified_refflat[["chrom", "exonStarts", "exonEnds", "name", "score", "strand", "exontype", "exon_position"]]
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