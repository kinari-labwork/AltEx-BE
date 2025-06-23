import pandas as pd
import pybedtools

def annotate_sequence_to_bed(bed: pd.DataFrame,species: str) -> pd.DataFrame:
    """
    Purpose:
        BED形式のデータに指定される遺伝子座位を参照して、指定された生物種の塩基配列をFASTAから取得し、dataframeに追加する
    Parameters:
        bed: pybedtools.Bedtools, BED形式のデータ
        species: str, 対象の生物種
        data: pd.DataFrame, アノテーションを与えるデータフレーム
    Returns:
        pd.DataFrame, アノテーションされたデータ
    """