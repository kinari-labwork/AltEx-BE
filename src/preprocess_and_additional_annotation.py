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
    import re
    # 正規表現パターンを使用して、染色体名が数字またはX, Yで終わるものを抜き出す（_random,_alt,_fixは除外）
    pattern = re.compile(r"^chr(\d+|X|Y)$")
    data_filtered = data[data["chrom"].str.match(pattern)]
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

def add_exon_position_flags(data: pd.DataFrame)-> pd.DataFrame:
    """
    Purpose:
        exon_position列を作成し、各行の転写産物に対してエキソンの位置を付与する
        各エキソンに'first','internal','last'のカテゴリを付加する
        エキソンが一つの場合は'single'を付加する
        のちにSA/SDを編集するsgRNAを作成するとき、1番目のエキソンのSA、最後のエキソンのSDを編集する意味がないから、事前にflagをつけておく
    Parameters:
        data: pd.DataFrame, refflatのデータフレーム
    """
    # 位置に応じて値を付与する関数の作成
    def get_category_list(x):
        n = len(x)
        if n ==1:
            return ['single']
        else:
            return ['first'] + ['internal'] * (n - 2) + ['last']
    data=data.copy()
    data["exon_position"] = data["exonStarts"].apply(get_category_list)
    return data