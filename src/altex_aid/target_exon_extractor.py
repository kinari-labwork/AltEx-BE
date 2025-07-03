from __future__ import annotations

import pandas as pd

import uuid

# BED形式も0base-start, 1base-endであるため、refFlatのexonStartsとexonEndsをそのまま使用する


def extract_target_exon(data: pd.DataFrame) -> pd.DataFrame:
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
    data = data[
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
    data = data.explode(["exonStarts", "exonEnds", "exontype", "exon_position"])
    data[["exonStarts", "exonEnds"]] = data[["exonStarts", "exonEnds"]].astype(
        int
    )  # int型に変換
    # exontypeがskippedまたはuniqueのエキソンだけを抽出
    data = data[data["exontype"].apply(lambda x: "skipped" in x or "unique" in x)]
    # 重複を削除し一方だけ残す
    data = data.drop_duplicates(subset=["exonStarts", "exonEnds"])
    data['name'] = [uuid.uuid4().hex for _ in range(len(data))]  # 一意のIDを生成
    data['score'] = 0  # BED形式のスコア列を追加
    #BED に合わせたカラム順に並べ替え
    data = data[["chrom", "exonStarts", "exonEnds", "name", "score", "strand", "exontype", "exon_position"]]
    return data.reset_index(drop=True)


def extract_splice_acceptor_regions(data: pd.DataFrame, window: int) -> pd.DataFrame:
    """
    Purpose :
        抜き出したskipped or unique exonのexonStart/Endから、SA部位周辺の、windowで指定した幅の座位を示すDataFrameを作成する
        strandが+の時はexonStartがSplice Acceptor, -の時はその逆でexonEndがSAになる
    """
    sa_data = data
    sa_data["chromStart"] = sa_data.apply(
        lambda row: row["exonStarts"] - window
        if row["strand"] == "+"
        else row["exonEnds"] - window,
        axis=1,
    )
    sa_data["chromEnd"] = sa_data.apply(
        lambda row: row["exonStarts"] + window
        if row["strand"] == "+"
        else row["exonEnds"] + window,
        axis=1,
    )
    return sa_data[["chrom","chromStart","chromEnd","name","score","strand"]].reset_index(drop=True)


def extract_splice_donor_regions(data: pd.DataFrame, window: int) -> pd.DataFrame:
    """
    Purpose :
        抜き出したskipped or unique exonのexonStart/Endから、SD部位周辺の、windowで指定した幅の座位を示すDataFrameを作成する
        strandが+の時はexonEndがSplice Donor, -の時はその逆でexonStartがSDになる
    """
    sd_data = data
    sd_data["chromStart"] = sd_data.apply(
        lambda row: row["exonEnds"] - window
        if row["strand"] == "+"
        else row["exonStarts"] - window,
        axis=1,
    )
    sd_data["chromEnd"] = sd_data.apply(
        lambda row: row["exonEnds"] + window
        if row["strand"] == "+"
        else row["exonStarts"] + window,
        axis=1,
    )
    return sd_data[["chrom","chromStart","chromEnd","name","score","strand"]].reset_index(drop=True)