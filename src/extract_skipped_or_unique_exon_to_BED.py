from __future__ import annotations
import pandas as pd

# BED形式も0base-start, 1base-endであるため、refFlatのexonStartsとexonEndsをそのまま使用する


def extract_skipped_or_unique_exon(data: pd.DataFrame)-> pd.DataFrame:
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
    data = data[["geneName","name", "chrom", "strand", "exonStarts", "exonEnds", "exontype"]].copy() 
    # 編集のために、リストになっている列を展開する
    data = data.explode(["exonStarts",'exonEnds', "exontype"])
    data[["exonStarts","exonEnds"]] = data[["exonStarts","exonEnds"]].astype(int) # int型に変換
    # exontypeがskippedまたはuniqueのエキソンだけを抽出
    data = data[data["exontype"].apply(lambda x: "skipped" in x or "unique" in x)]
    # 重複を削除し一方だけ残す
    data = data.drop_duplicates(subset=["exonStarts", "exonEnds"])
    # index列を追加
    data = data.reset_index(drop=True)
    data["index"] = data.index
    return data.reset_index(drop=True)

def format_to_single_exon_bed(data: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        1エキソン1行のデータフレームをBED形式に変換する
    Parameters:
        data: pd.DataFrame, 1エキソン1行のデータフレーム
    Returns:
        pd.DataFrame, BED形式に変換されたデータフレーム (BEDはタブ区切り形式なので実行時にタブ区切りに変更して保存する必要がある)
    """
    bed_data = data[["chrom", "exonStarts", "exonEnds","geneName","strand","index"]].copy()
    

    bed_data.columns = ["chrom", "chromStart", "chromEnd", "name", "strand","score"]
    return bed_data.reset_index(drop=True)
