import pandas as pd
import mappy as mp
import pathlib as Path
from altex_aid.sgrna_designer import BaseEditor


def convert_empty_list_into_na(target_exon_with_sgrna_dict: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Convert empty lists in the DataFrame to NaN.
    """
    for df in target_exon_with_sgrna_dict.values():
        for col in df.columns:
            df[col] = df[col].apply(lambda x: pd.NA if isinstance(x, list) and not x else x)
    return df

def prepare_melted_df(target_exon_with_sgrna_df: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose : 1行-複数のsgRNAの状態からexplodeして1行-1sgRNAに変換する前処理として、列名を整形し、データフレームを整形する
    Parameters : target_exon_with_sgrna_df: あるBEに対して設計されたsgRNAのdf。 acceptorとdonorの2つのサイトタイプの情報が含まれる
    Return : melted_df: 整形後のデータフレーム
    """

    foundation_cols = ["chrom", "exonStarts", "exonEnds", "strand", "exontype", "exon_position", "base_editor_name"]
    df_list = []
    for site in ["acceptor", "donor"]:
        cols = [col for col in target_exon_with_sgrna_df.columns if col.startswith(site)]
        site_df = target_exon_with_sgrna_df[foundation_cols + cols].rename(
            columns={
                col: f"{col.replace(f'{site}_', '')}" for col in cols
            }
        )
        site_df = site_df.dropna(subset=["sgrna_target_sequence"])
        site_df["site_type"] = site
        df_list.append(site_df)
    melted_df = pd.concat(df_list, ignore_index=True)
    return melted_df

def validate_exploded_df(melted_df: pd.DataFrame) -> None:
    """
    Purpose: Validate the exploded sgRNA DataFrame.
    """
    if melted_df.empty:
        raise ValueError("there are no designed sgRNAs for your interest gene and your base editor.")

def explode_sgrna_df(target_exon_with_sgrna_dict: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Purpose: 複数のBEに対して設計されたsgRNAのdfが、1列-1sgRNAになるように変換する
    """
    # 各BEに対して、sgRNAのdfを1列-1sgRNAに変換
    exploded_dfs = []
    for be, df in target_exon_with_sgrna_dict.items():
        df["base_editor_name"] = be
        df = prepare_melted_df(df)
        sgrna_cols = [col for col in df.columns if col.startswith("sgrna_")]
        df = df.explode(column=sgrna_cols)
        exploded_dfs.append(df)
    # すべてのBEのdfを結合
    exploded_sgrna_df = pd.concat(exploded_dfs, ignore_index=True)
    return exploded_sgrna_df

def add_base_editor_info_to_row(row: pd.Series, base_editors: list[BaseEditor]) -> pd.Series:
    """
    Purpose: Add base editor information to the exploded sgRNA DataFrame.
    """
    for be in base_editors:
        if row["base_editor_name"] == be.base_editor_name:
            row["base_editor_pam"] = be.pam_sequence
            row["base_editor_editing_window_start"] = be.editing_window_start_in_grna
            row["base_editor_editing_window_end"] = be.editing_window_end_in_grna
            row["base_editor_type"] = be.base_editor_type
            break
    return row

def add_base_editor_info_to_df(exploded_sgrna_df: pd.DataFrame, base_editors: list[BaseEditor]) -> pd.DataFrame:
    """
    Purpose: Add base editor information to the exploded sgRNA DataFrame.
    """
    base_editor_name_col_index = exploded_sgrna_df.get_loc("base_editor_name") + 1
    exploded_sgrna_df.insert(base_editor_name_col_index, ["base_editor_pam", "base_editor_editing_window_start", "base_editor_editing_window_end", "base_editor_type"], pd.NA)
    exploded_sgrna_df = exploded_sgrna_df.apply(lambda row: add_base_editor_info_to_row(row, base_editors), axis=1)
    return exploded_sgrna_df

def add_crisprdirect_url_to_row(row: pd.Series, assembly_name:str) -> pd.Series:
    """
    Purpose : 列ごとにCRISPRdirectのURLを追加する
    """
    target_sequence = row["sgrna_target_sequence"].replace('+', '').lower()
    base_url = "https://crispr.dbcls.jp/?userseq="
    row["crisprdirect_url"] = base_url + target_sequence + "&pam=" + row["base_editor_pam"] + "&db=" + assembly_name
    return row

def add_crisprdirect_url_to_df(exploded_sgrna_df: pd.DataFrame, assembly_name: str) -> pd.DataFrame:
    """
    Purpose : DataFrameの各行にCRISPRdirectのURLを追加する
    """
    exploded_sgrna_df = exploded_sgrna_df.apply(lambda row: add_crisprdirect_url_to_row(row, assembly_name), axis=1)
    return exploded_sgrna_df

def calculate_offtarget_site_count(row:pd.Series, aligner:mp.Aligner, fasta_path:Path) -> pd.Series:
    """
    Purpose: CRISPR direct に飛ぶ前に簡易的にPAM+20bpの完全一致のオフターゲット数を計算する
    """
    aligner = mp.Aligner(str(fasta_path))
    pam_plus_target_seq = row["pam_plus_target_sequence"].replace('+', '').lower()

    exact_match_count = 0
    for name, seq in mp.fastx_read(fasta_path):
        for hit in aligner.map(pam_plus_target_seq):
            if hit.mlen == len(pam_plus_target_seq):
                exact_match_count += 1
            if exact_match_count > 10: # 検索時間を抑えるため10個以上のマッチが見つかったら、計算を中止
                break

    row["exact_match_count_of_20bp_plus_pam"] = exact_match_count

    return row
