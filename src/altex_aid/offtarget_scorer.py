import pandas as pd
import mappy as mp
import pathlib as Path

def add_crisprdirect_url_to_df(exploded_sgrna_df: pd.DataFrame, assembly_name: str) -> pd.DataFrame:
    """
    Purpose: exploded_sgrna_dfにCRISPRdirectのURLを追加する
    """
    base_url = "https://crispr.dbcls.jp/?userseq="
    
    # 列全体に対して一度に文字列操作を行う
    target_sequences = exploded_sgrna_df["sgrna_target_sequence"].str.replace('+', '', regex=False).str.lower()
    pams = exploded_sgrna_df["base_editor_pam"] # 事前にマージしておく必要がある
    
    exploded_sgrna_df["crisprdirect_url"] = base_url + target_sequences + "&pam=" + pams + "&db=" + assembly_name
    return exploded_sgrna_df


def calculate_offtarget_site_count_to_df(exploded_sgrna_df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
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
