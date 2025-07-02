from __future__ import annotations

import pandas as pd

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
    data = data[data["exontype"].apply(lambda x: "skipped" in x or "unique" in x)]
    data = data[
        [
            "geneName",
            "name",
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
    # index=(score)列を追加 (BEDとして扱いやすくするため、名称をscoreにしておく)
    data = data.reset_index(drop=True)
    data["score"] = data.index
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
    return sa_data[
        ["geneName", "chrom", "strand", "chromStart", "chromEnd", "score"]
    ].reset_index(drop=True)


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
    return sd_data[
        ["geneName", "chrom", "strand", "chromStart", "chromEnd", "score"]
    ].reset_index(drop=True)


def format_to_single_exon_bed(data: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        1エキソン1行のデータフレームをBED形式に変換する
    Parameters:
        data: pd.DataFrame, 1エキソン1行のデータフレーム
    Returns:
        pd.DataFrame, BED形式に変換されたデータフレーム (BEDはタブ区切り形式なので実行時にタブ区切りに変更して保存する必要がある)
    """
    bed_data = data[["chrom", "chromStart", "chromEnd", "geneName", "score", "strand"]]
    bed_data["name"] = bed_data["geneName"]  # geneNameをnameとして使用
    # BED形式を満たすように列を並べ替える
    bed_data = bed_data[["chrom", "chromStart", "chromEnd", "name", "score", "strand"]]
    return bed_data.reset_index(drop=True)
