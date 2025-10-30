from __future__ import annotations

import re
import pandas as pd
import logging
from . import logging_config # noqa: F401

def select_interest_genes(refFlat: pd.DataFrame, interest_genes: list[str]) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームから、興味のある遺伝子のみを選択する。
    Parameters:
        refFlat: pd.DataFrame, refFlatのデータフレーム
        interest_genes: list[str], 興味のある遺伝子名のリスト(gene symbol または Refseq ID)
    Returns:
        pd.DataFrame, 興味のある遺伝子のみを含むrefFlatのデータフレーム
    """
    gene_symbol_set = set(refFlat["geneName"].values)
    ref_seq_id_set = set(refFlat["name"].values)

    for gene in interest_genes:
        if gene not in gene_symbol_set and gene not in ref_seq_id_set:
            logging.warning(f"Gene {gene} is not found in refFlat.")
            continue
        else :
            logging.info(f"Gene {gene} is found in refFlat.")
    
    refFlat = refFlat[refFlat["geneName"].isin(interest_genes) | refFlat["name"].isin(interest_genes)].reset_index(drop=True)
    return refFlat

def check_multiple_exon_existance(refFlat: pd.DataFrame, interest_gene_list) -> bool:
    """
    Purpose:
        refFlatのデータフレームに、複数のエキソンが存在するかを確認する。
    Parameters:
        refFlat: pd.DataFrame, refFlatのデータフレーム
    Returns:
        bool, 複数のエキソンが存在する場合はTrue、存在しない場合はFalse
    """
    found = False
    for gene in interest_gene_list:
        exon_counts = refFlat[refFlat["geneName"] == gene]["exonCount"]
        if (exon_counts > 1).any():
            logging.info(f"Gene {gene} has multiple exons")
            found = True
    return found

# constitutive exonも含めてデザインするなら、いらない可能性もある
def check_transcript_variant(refFlat: pd.DataFrame, interest_genes: list[str]) -> bool:
    """
    Purpose:
        refFlatのデータフレームに、トランスクリプトのバリアントが存在するかを確認する。
    Parameters:
        refFlat: pd.DataFrame, refFlatのデータフレーム
    Returns:
        bool, トランスクリプトのバリアントが存在する場合はTrue、存在しない場合はFalse
    """
    bool_list = []
    for gene in interest_genes:
        # 遺伝子ごとにトランスクリプトの数をカウント
        transcripts = refFlat[refFlat["geneName"] == gene]
        if transcripts.shape[0] > 1:
            logging.info(f"Gene {gene} has multiple transcripts")
            bool_list.append(True)
        else:
            logging.warning(f"Gene {gene} has a single transcript")
            bool_list.append(False)
    if all([x is False for x in bool_list]):
        logging.warning("All genes have a single transcript, stop further processing.")
        return False
    return True


def parse_exon_coordinates(refFlat: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        exonStart exonEndが別々のカラムに格納されているので、(start, end)のタプルのリストに変換する。
    Parameters:
        refFlat: pd.DataFrame, refFlatのデータフレーム
    """
    # Convert the exonStarts and exonEnds columns to lists of integers
    refFlat["exonStarts"] = refFlat["exonStarts"].apply(
        lambda x: [int(i) for i in x.split(",") if i.strip() != ""]
    )
    refFlat["exonEnds"] = refFlat["exonEnds"].apply(
        lambda x: [int(i) for i in x.split(",") if i.strip() != ""]
    )

    refFlat["exons"] = refFlat.apply(
        lambda row: list(zip(row["exonStarts"], row["exonEnds"])), axis=1
    )
    return refFlat


def calculate_exon_lengths(refFlat: pd.DataFrame) -> pd.DataFrame:
    """Purpose:
        refFlatのデータフレームに、各エキソンの長さを計算して追加する。
    Parameters:
        refFlat: pd.DataFrame, refFlatのデータフレーム
    Returns:
        pd.DataFrame, 各エキソンの長さを追加したrefFlatのデータフレーム
    """
    # Calculate the lengths of each exon
    # refflatのstartは0-baseでendは1-baseなので、毎回1を足す必要がない
    refFlat["exonlengths"] = refFlat.apply(
        lambda row: [
            end - start for start, end in zip(row["exonStarts"], row["exonEnds"])
        ],
        axis=1,
    )
    return refFlat


def drop_abnormal_mapped_transcripts(refflat: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームから、異常な染色体にマッピングされたトランスクリプトを削除する。
    Parameters:
        refflat: pd.DataFrame, refFlatのデータフレーム
    Returns:
        pd.DataFrame, 異常な染色体マッピングを持つトランスクリプトを削除したrefFlatのデータフレーム
    """

    # 正規表現パターンを使用して、染色体名が数字またはX, Yで終わるものを抜き出す（_random,_alt,_fixは除外）
    pattern = re.compile(r"^chr(\d+|X|Y)$")
    data_filtered = refflat[refflat["chrom"].str.match(pattern)]
    return data_filtered.reset_index(drop=True)


def annotate_cording_information(refflat: pd.DataFrame) -> pd.DataFrame:
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
    refflat["coding"] = ""
    refflat.loc[refflat["name"].str.match(cording_pattern), "coding"] = "coding"
    refflat.loc[refflat["name"].str.match(non_coding_pattern), "coding"] = "non-coding"
    refflat["coding"] = refflat["coding"].astype("category")
    return refflat


def annotate_flame_information(refflat: pd.DataFrame) -> pd.DataFrame:
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
        return ["in-flame" if length % 3 == 0 else "out-flame" for length in lengths]

    refflat["flame"] = refflat["exonlengths"].apply(calc_flame)
    return refflat


def annotate_variant_count(refflat: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームに、バリアント数を追加する。
    Parameters:
        data: pd.DataFrame, refFlatのデータフレーム
        variant_data: pd.DataFrame, バリアント情報のデータフレーム
    Returns:
        pd.DataFrame, バリアント数を追加したrefFlatのデータフレーム
    """
    variant_counts = refflat.groupby("geneName")["name"].nunique().reset_index()
    variant_counts.columns = ["geneName", "variant_count"]
    refflat = refflat.merge(variant_counts, on="geneName", how="left")
    return refflat


def add_exon_position_flags(refflat: pd.DataFrame) -> pd.DataFrame:
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
        if n == 1:
            return ["single"]
        else:
            return ["first"] + ["internal"] * (n - 2) + ["last"]

    refflat["exon_position"] = refflat["exonStarts"].apply(get_category_list)

    def flip_first_last_to_minus_strand(row):
        """
        Purpose:
        マイナス鎖の転写産物に対して、firstとlastを入れ替える
        なぜなら、マイナス鎖の転写産物では、最初のエキソンが最後のエキソンになり、最後のエキソンが最初のエキソンになるから
        Parameters:
            row: pd.Series, 各行のデータ
        """
        if row["strand"] == "-":
            row["exon_position"] = row["exon_position"][::-1]
        return row

    refflat = refflat.apply(flip_first_last_to_minus_strand, axis=1)
    return refflat

def annotate_utr_and_cds_exons(refflat: pd.DataFrame) -> pd.DataFrame:
    """
    各エキソンごとに 'cds_exon', 'cds_edge_exon', 'utr_exon' のラベルを付与し、cds_infoカラムに格納する。
    """
    def label_exons(row):
        # non-coding遺伝子はすべてutr_exon
        if row["coding"] == "non-coding":
            return ["utr_exon" for _ in row["exons"]]
        cds_start = row["cdsStart"]
        cds_end = row["cdsEnd"]
        exon_starts = row["exonStarts"]
        exon_ends = row["exonEnds"]
        labels = []
        for start, end in zip(exon_starts, exon_ends):
            if end <= cds_start or start >= cds_end:
                label = "utr_exon"
            elif (start < cds_start < end) and (start < cds_end < end):
                label = "cds_edge_exon_start_end"
            elif start < cds_start < end:
                label = "cds_edge_exon_start"
            elif start < cds_end < end:
                label = "cds_edge_exon_end"
            else:
                label = "cds_exon"
            labels.append(label)
        return labels

    refflat["cds_info"] = refflat.apply(label_exons, axis=1)
    return refflat

def preprocess_refflat(refflat: pd.DataFrame, interest_genes: list[str]) -> pd.DataFrame:
    """
    このモジュールの関数をwrapした関数
    """
    refflat = select_interest_genes(refflat, interest_genes)
    if refflat.empty:
        return refflat
    refflat = parse_exon_coordinates(refflat)
    refflat = calculate_exon_lengths(refflat)
    refflat = drop_abnormal_mapped_transcripts(refflat)
    refflat = annotate_cording_information(refflat)
    refflat = annotate_flame_information(refflat)
    refflat = add_exon_position_flags(refflat)
    refflat = annotate_utr_and_cds_exons(refflat)

    return refflat