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
        melted_df = prepare_melted_df(df)
        sgrna_cols = [col for col in melted_df.columns if col.startswith("sgrna_")]
        melted_df = melted_df.explode(column=sgrna_cols)
        exploded_dfs.append(melted_df)
    # すべてのBEのdfを結合
    exploded_sgrna_df = pd.concat(exploded_dfs, ignore_index=True)
    return exploded_sgrna_df

def add_base_editor_info_to_df_revised(exploded_sgrna_df: pd.DataFrame, base_editors: list[BaseEditor]) -> pd.DataFrame:
    """
    Puropose: exploded_sgrna_dfにBaseEditorの情報を追加する
    """
    # BaseEditorのリストをDataFrameに変換
    be_info = [
        {
            "base_editor_name": be.base_editor_name,
            "base_editor_pam": be.pam_sequence,
            "base_editor_editing_window_start": be.editing_window_start_in_grna,
            "base_editor_editing_window_end": be.editing_window_end_in_grna,
            "base_editor_type": be.base_editor_type,
        }
        for be in base_editors
    ]
    be_df = pd.DataFrame(be_info)

    # 'base_editor_name'をキーとしてマージ（結合）
    exploded_sgrna_df = pd.merge(exploded_sgrna_df, be_df, on="base_editor_name", how="left")
    return exploded_sgrna_df

def add_crisprdirect_url_to_df_revised(exploded_sgrna_df: pd.DataFrame, assembly_name: str) -> pd.DataFrame:
    """
    Purpose: exploded_sgrna_dfにCRISPRdirectのURLを追加する
    """
    base_url = "https://crispr.dbcls.jp/?userseq="
    
    # 列全体に対して一度に文字列操作を行う
    target_sequences = exploded_sgrna_df["sgrna_target_sequence"].str.replace('+', '', regex=False).str.lower()
    pams = exploded_sgrna_df["base_editor_pam"] # 事前にマージしておく必要がある
    
    exploded_sgrna_df["crisprdirect_url"] = base_url + target_sequences + "&pam=" + pams + "&db=" + assembly_name
    return exploded_sgrna_df


def calculate_offtarget_site_count_to_df_revised(exploded_sgrna_df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    """
    Purpose: exploded_sgrna_dfにオフターゲットの数を追加する
    """
    aligner = mp.Aligner(str(fasta_path))
    
    # 計算対象のユニークな配列を取得
    unique_sequences = exploded_sgrna_df["sgrna_target_sequence"].dropna().unique()

    # ユニークな配列に対してオフターゲット数を計算し、辞書に保存
    offtarget_counts = {}
    for seq in unique_sequences:
        # mappyの正しい使い方。内部でインデックスが使われるため、FASTAを再読み込みする必要はない。
        # また、"+"は塩基ではないため、シーケンスから除外すべき
        clean_seq = seq.replace('+', '', 1).lower()
        count = sum(1 for hit in aligner.map(clean_seq) if hit.mlen == len(clean_seq))
        offtarget_counts[seq] = count

    # 計算結果を元のDataFrameにマップする
    exploded_sgrna_df["exact_match_count_of_20bp_plus_pam"] = exploded_sgrna_df["pam_plus_target_sequence"].map(offtarget_counts)
    return exploded_sgrna_df
