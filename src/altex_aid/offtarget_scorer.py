import pandas as pd
from pathlib import Path
import re
import logging
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

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

def calculate_offtarget_site_count_simple(exploded_sgrna_df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    """
    Purpose: exploded_sgrna_dfにオフターゲットサイトのカウントを追加する
    Parameters: exploded_sgrna_df (pd.DataFrame): 入力のDataFrame
                fasta_path (Path): FASTAファイルのパス
    Return: pd.DataFrame: オフターゲットサイトのカウントを追加したDataFrame
    """
    # + を除去し、全て大文字に変換
    sgrna_sequences = exploded_sgrna_df["sgrna_target_sequence"].str.replace('+', '', regex=False).str.upper()
    unique_sequences = sgrna_sequences.dropna().unique()
    logging.info(f"count off-target site of {len(unique_sequences)} sgRNAs")
    compiled_patterns = {seq: re.compile(f"(?={seq})") for seq in unique_sequences}
    offtarget_count_dict = {seq:0 for seq in unique_sequences}
    with open(fasta_path, 'r') as fasta_file:
        chrom_count = sum(1 for line in fasta_file if line.startswith(">"))
    logging.info(f"detected number of chromosome is : {chrom_count}")

    with open(fasta_path, 'r') as fasta_file:
        current_chromosome = []
        chromosome_sequence = ""
        chrom_name = None
        pbar_chrom = tqdm(total=chrom_count, desc="Chromosome", position=0)
        for line in fasta_file:
            if line.startswith(">"):
                if current_chromosome:
                    chromosome_sequence = "".join(current_chromosome)
                    # sgRNAごとの進捗バー
                    for seq in tqdm(unique_sequences, desc=f"sgRNA in {chrom_name}", position=1, leave=False):
                        pattern = compiled_patterns[seq]
                        offtarget_count_dict[seq] += len(list(pattern.finditer(chromosome_sequence)))
                current_chromosome = []
                chrom_name = line[1:].strip()
                pbar_chrom.update(1)
            else:
                # FASTA側も大文字に変換する
                current_chromosome.append(line.strip().upper())
        if current_chromosome:
            for seq in tqdm(unique_sequences, desc=f"sgRNA in {chrom_name}", position=1, leave=False):
                pattern = compiled_patterns[seq]
                offtarget_count_dict[seq] += len(list(pattern.finditer("".join(current_chromosome))))
        pbar_chrom.close()

    exploded_sgrna_df["pam+20bp_exact_match_count"] = sgrna_sequences.map(offtarget_count_dict)
    return exploded_sgrna_df

def score_offtargets(exploded_sgrna_df: pd.DataFrame, assembly_name: str, fasta_path: Path) -> pd.DataFrame:
    """
    Purpose: このモジュールのラップ関数
    """
    exploded_sgrna_df = add_crisprdirect_url_to_df(exploded_sgrna_df, assembly_name)
    exploded_sgrna_df = calculate_offtarget_site_count_simple(exploded_sgrna_df, fasta_path)
    return exploded_sgrna_df