import pandas as pd
import mappy as mp
from pathlib import Path

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

def calculate_offtarget_site_count(exploded_sgrna_df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    """
    Purpose: とりあえず 20bp + PAMの完全一致のマッチ数を計算する
    """
    aligner = mp.Aligner(str(fasta_path), preset="sr")

    # 1. 計算対象の列を準備（+を削除し、小文字に統一）
    sgrna_sequences = exploded_sgrna_df["sgrna_target_sequence"].str.replace('+', '', regex=False).str.lower()

    unique_sequences = sgrna_sequences.dropna().unique()
    # 2. ユニークな各配列に対してオフターゲット数を計算し、結果を辞書に保存
    offtarget_counts = {}
    for seq in unique_sequences:
        exact_match_count = 0
        for hit in aligner.map(seq):
            if hit.mlen == len(seq) and hit.NM == 0:
                exact_match_count += 1
            if exact_match_count > 10:
                break
        offtarget_counts[seq] = exact_match_count

    # 3. 計算結果を元のDataFrameにマップする
    exploded_sgrna_df["pam+20bp_exact_match_count"] = sgrna_sequences.map(offtarget_counts)

    return exploded_sgrna_df
