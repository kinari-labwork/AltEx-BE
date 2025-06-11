from __future__ import annotations
import pandas as pd

def drop_abnormal_mapped_transcripts(data: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームから、異常な染色体にマッピングされたトランスクリプトを削除する。
    Parameters:
        refflat: pd.DataFrame, refFlatのデータフレーム
    Returns:
        pd.DataFrame, 異常な染色体マッピングを持つトランスクリプトを削除したrefFlatのデータフレーム
    """
    mask = data["chrom"].str.endswith("_random") | data["chrom"].str.startswith("chrUn") | data["chrom"].str.endswith("_alt")
    data_filtered = data[~mask]  # ~mask で「条件に一致しない」行を選ぶ
    return data_filtered.reset_index(drop=True)

def cording_information_annotator(data: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームに、コーディング情報を追加する。
    Parameters:
        data: pd.DataFrame, refFlatのデータフレーム
    Returns:
        pd.DataFrame, コーディング情報を追加したrefFlatのデータフレーム
    """
    # コーディングと非コーディングのトランスクリプトを識別するための正規表現パターン
    # NMはコーディング、NRは非コーディング
    import re
    cording_pattern = re.compile(r"^NM")
    non_coding_pattern = re.compile(r"^NR")
    data["coding"] ="" 
    data.loc[data["name"].str.match(cording_pattern), "coding"] = "coding"
    data.loc[data["name"].str.match(non_coding_pattern), "coding"] = "non-coding"
    data["coding"] = data["coding"].astype("category")
    return data

def flame_information_annotator(data: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームに、フレーム情報を追加する。
    Parameters:
        data: pd.DataFrame, refFlatのデータフレーム
    Returns:
        pd.DataFrame, フレーム情報を追加したrefFlatのデータフレーム
    """
    # exonlengths列（リスト）に対してmod3を計算し、0ならin-flame, それ以外はout-flame
    def calc_flame(lengths):
        return ["in-flame" if l % 3 == 0 else "out-flame" for l in lengths]

    data = data.copy()
    data["flame"] = data["exonlengths"].apply(calc_flame)
    return data


def variant_count_annotator(data: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームに、バリアント数を追加する。
    Parameters:
        data: pd.DataFrame, refFlatのデータフレーム
        variant_data: pd.DataFrame, バリアント情報のデータフレーム
    Returns:
        pd.DataFrame, バリアント数を追加したrefFlatのデータフレーム
    """
    variant_counts = data.groupby("geneName")["name"].nunique().reset_index()
    variant_counts.columns = ["geneName", "variant_count"]
    data = data.merge(variant_counts, on="geneName", how="left")
    return data