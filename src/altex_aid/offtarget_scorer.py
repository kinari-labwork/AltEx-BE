import pandas as pd
from pathlib import Path
import logging
from tqdm import tqdm
import ahocorasick
from . import logging_config # noqa: F401

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

def convert_dna_to_reversed_complement_dna(sequence: str) -> str:
    """
    purpose:
        塩基配列を逆相補のRNAに変換する
    Parameters:
        sequence: 変換したいDNA配列
    Returns:
        reversed_complement_rna_sequence: 入力したDNA配列の逆相補RNA配列
    """
    complement_map = {
        "A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
        "a": "t", "t": "a", "c": "g", "g": "c", "+": "+"
    }
    return "".join([complement_map[base] for base in reversed(sequence)])

def add_reversed_complement_sgrna_column(exploded_sgrna_df: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose: 逆相補のsgRNA配列を追加する
    """
    exploded_sgrna_df["reversed_sgrna_target_sequence"] = exploded_sgrna_df["sgrna_target_sequence"].apply(
        convert_dna_to_reversed_complement_dna
    )
    return exploded_sgrna_df

def calculate_offtarget_site_count_ahocorasick(exploded_sgrna_df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    # 遺伝子が - strandの場合、出力されている配列は - strandの配列である。しかし、検索対象は+ strandであるため、逆相補に変換する必要がある。
    unique_sequences = set()
    unique_sequences = exploded_sgrna_df.apply(
        lambda row: row["reversed_sgrna_target_sequence"].replace('+', '').upper() if row["strand"] == '-' else row["sgrna_target_sequence"].replace('+', '').upper(),
        axis=1
    )

    # まず最初にAho-CorasickのAutomatonを構築
    automaton = ahocorasick.Automaton()
    for idx, seq in enumerate(unique_sequences):
        automaton.add_word(seq, (idx, seq))
    automaton.make_automaton()

    # sgRNAごとのカウント辞書
    offtarget_count_dict = {seq: 0 for seq in unique_sequences}

    with open(fasta_path, 'r') as fasta_file:
        header_count = sum(1 for line in fasta_file if line.startswith(">"))
        logging.info(f"Number of chromosomes in your FASTA file: {header_count}")
        pbar = tqdm(total=header_count, desc="Calculating off-target counts", unit="chromosome")
        fasta_file.seek(0)
        chrom_seq = ""
        for line in fasta_file:
            if line.startswith(">"):
                # 新しい染色体に切り替え
                if chrom_seq:
                    chrom_seq = chrom_seq.upper()
                    for end_idx, (idx, seq) in automaton.iter(chrom_seq):
                        offtarget_count_dict[seq] += 1
                    chrom_seq = ""
                    pbar.update(1)
            else:
                chrom_seq += line.strip()
        # 最後の染色体も処理
        if chrom_seq:
            chrom_seq = chrom_seq.upper()
            for end_idx, (idx, seq) in automaton.iter(chrom_seq):
                offtarget_count_dict[seq] += 1
            pbar.update(1)
            pbar.close()

    exploded_sgrna_df["pam+20bp_exact_match_count"] = exploded_sgrna_df.apply(
        lambda row: offtarget_count_dict.get(
            (row["reversed_sgrna_target_sequence"] if row["strand"] == '-' else row["sgrna_target_sequence"])
            .replace('+', '').upper(),
            0  # デフォルト値
    ),
    axis=1
    )
    exploded_sgrna_df = exploded_sgrna_df.drop(columns=["reversed_sgrna_target_sequence"])
    return exploded_sgrna_df

def score_offtargets(exploded_sgrna_df: pd.DataFrame, assembly_name: str, fasta_path: Path) -> pd.DataFrame:
    """
    Purpose: このモジュールのラップ関数
    """
    exploded_sgrna_df = add_crisprdirect_url_to_df(exploded_sgrna_df, assembly_name)
    exploded_sgrna_df = add_reversed_complement_sgrna_column(exploded_sgrna_df)
    exploded_sgrna_df = calculate_offtarget_site_count_ahocorasick(exploded_sgrna_df, fasta_path)
    return exploded_sgrna_df