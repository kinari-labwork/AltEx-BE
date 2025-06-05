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
    mask = data["chrom"].str.startswith("chrX_GL") | data["chrom"].str.startswith("chr1_GL") | data["chrom"].str.startswith("chr1_MU") | data["chrom"].str.startswith("chr1_Un")
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
    data["coding"] ="" 
    # cdsStartとcdsEndが同じ場合はnon-coding、異なる場合はcodingとする
    data.loc[data["cdsStart"] == data["cdsEnd"], "coding"] = "non-coding" 
    data.loc[data["cdsStart"] != data["cdsEnd"], "coding"] = "coding"
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
    data["flame"] = data["exonlength"].apply(calc_flame)
    return data

def max_min_exon_count_annotator(data: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームに、最大および最小のエクソン数を追加する。
    Parameters:
        data: pd.DataFrame, refFlatのデータフレーム
    Returns:
        pd.DataFrame, 最大および最小のエクソン数を追加したrefFlatのデータフレーム
    """
    maximum_exon_counts = data.groupby("geneName")["exonCount"].max().reset_index()
    maximum_exon_counts.columns = ["geneName", "max_exon_count"]
    minimum_exon_counts = data.groupby("geneName")["exonCount"].min().reset_index()
    minimum_exon_counts.columns = ["geneName", "min_exon_count"]
    data = data.merge(maximum_exon_counts, on="geneName", how="left")
    data = data.merge(minimum_exon_counts, on="geneName", how="left")
    data["max_exon_count"] = data["max_exon_count"].astype("Int64")
    data["min_exon_count"] = data["min_exon_count"].astype("Int64")
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
    data["variant_count"] = data["variant_count"].astype("Int64")
    return data