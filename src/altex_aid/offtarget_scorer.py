import pandas as pd
import mappy as mp
from pathlib import Path
import io 

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

def make_fasta_of_sgrnas(exploded_sgrna_df: pd.DataFrame) -> str:
    """
    Purpose: Create a FASTA format string from the exploded_sgrna_df.
    """
    fasta_str = ""
    exploded_sgrna_df["sgrna_target_sequence"] = exploded_sgrna_df["sgrna_target_sequence"].str.replace('+', '', regex=False).str.lower()
    for _, row in exploded_sgrna_df.iterrows():
        fasta_str += f"> {row['uuid']}\n{row['sgrna_target_sequence']}\n"
    return fasta_str

def calculate_offtarget_site_count_to_df(exploded_sgrna_df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    """
    Purpose: exploded_sgrna_dfにオフターゲットの数を追加する
    Parameters: exploded_sgrna_df: 入力のDataFrame, sgrna_fasta: sgRNAのFASTAファイル, fasta_path: 対象のFASTAファイルパス
    """
    aligner = mp.Aligner(str(fasta_path))
    exploded_sgrna_fasta = make_fasta_of_sgrnas(exploded_sgrna_df)

    match_dict = {}

    fasta_io = io.StringIO(exploded_sgrna_fasta)
    for name, seq in mp.fastx_read(fasta_io): 
        exact_match_count = 0
        for hit in aligner.map(seq):
            if hit.mlen == len(seq) and hit.NM == 0:
                exact_match_count += 1
            if exact_match_count > 10:
                match_dict[name] = exact_match_count
                break
        match_dict[name] = exact_match_count

    exploded_sgrna_df["pam+20bp_exact_match_count"] = exploded_sgrna_df["uuid"].map(match_dict).astype(int)
    return exploded_sgrna_df
