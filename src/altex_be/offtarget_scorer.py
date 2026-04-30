import pandas as pd
from pathlib import Path
import logging
from tqdm import tqdm
import ahocorasick
from itertools import product
from . import logging_config # noqa: F401

SEED_LENGTH = 12
IUPAC_TO_BASES = {
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
    "M": ("A", "C"),
    "R": ("A", "G"),
    "W": ("A", "T"),
    "S": ("C", "G"),
    "Y": ("C", "T"),
    "K": ("G", "T"),
    "V": ("A", "C", "G"),
    "H": ("A", "C", "T"),
    "D": ("A", "G", "T"),
    "B": ("C", "G", "T"),
    "N": ("A", "C", "G", "T"),
}

def add_crisprdirect_url_to_df(exploded_sgrna_df: pd.DataFrame, assembly_name: str) -> pd.DataFrame:
    """
    Purpose: exploded_sgrna_dfにCRISPRdirectのURLを追加する
    """
    base_url = "https://crispr.dbcls.jp/?userseq="
    
    # 列全体に対して一度に文字列操作を行う
    target_sequences = exploded_sgrna_df["sgrna_target_sequence"].str.replace('+', '', regex=False).str.lower()
    pams = exploded_sgrna_df["base_editor_pam_sequence"]
    
    exploded_sgrna_df["crisprdirect_url"] = base_url + target_sequences + "&pam=" + pams + "&db=" + assembly_name
    return exploded_sgrna_df

def convert_dna_to_reversed_complement_dna(sequence: str) -> str:
    """
    purpose:
        塩基配列を逆相補のDNAに変換する
    Parameters:
        sequence: 変換したいDNA配列
    Returns:
        reversed_complement_dna_sequence: 入力したDNA配列の逆相補DNA配列
    """
    complement_map = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
        "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
        "a": "t", "t": "a", "c": "g", "g": "c",
        "r": "y", "y": "r", "s": "s", "w": "w", "k": "m", "m": "k",
        "b": "v", "v": "b", "d": "h", "h": "d", "n": "n"
    }
    return "".join([complement_map[base] for base in reversed(sequence)])

def expand_pam_iupac(pam_rule: str) -> list[str]:
    """
    Purpose:
        IUPAC表記のPAM配列をすべての具体的な塩基配列に展開する。
    """
    pam_rule = pam_rule.upper()
    try:
        base_choices = [IUPAC_TO_BASES[base] for base in pam_rule]
    except KeyError as exc:
        raise ValueError(f"Unsupported IUPAC base in PAM rule: {pam_rule}") from exc
    return ["".join(bases) for bases in product(*base_choices)]

def get_seed_sequence_from_protospacer(protospacer: str, seed_len: int = SEED_LENGTH) -> str:
    """
    Purpose:
        protospacer配列（20bp想定）から、PAM隣接側seed領域を抽出する。
        protospacerの右側にPAMがある向きを基準とし、seedは3'側末端を使う。
    Parameters:
        protospacer: PAMなしの標的配列
        seed_len: Seed領域の長さ (default: 12)
    Returns:
        str: Seed領域配列
    """
    return protospacer[-seed_len:]

def calculate_offtarget_site_count_ahocorasick(exploded_sgrna_df: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    """
    Purpose : ahocorasick法を用いて PAM+20bpのオフターゲットサイト数を計算する
    Parameters : exploded_sgrna_df: sgRNAが1行1sgRNAに展開されたデータフレーム, fasta_path: FASTAファイルのパス
    Returns : exploded_sgrna_df: PAM+20bpのオフターゲットサイト数を追加したデータフレーム
    Algorism : sgRNA配列とその逆相補配列をセットに追加し、Aho-CorasickのAutomatonを構築。各染色体配列に対してAutomatonを用いて検索し、各sgRNA配列の出現回数をカウントする。
    """
    full_sequences_per_row = []
    seed_sequences_per_row = []
    all_full_sequences = set()
    all_seed_sequences = set()

    for _, row in exploded_sgrna_df.iterrows():
        protospacer = row["sgrna_sequence"].upper()
        pam_rule = row["base_editor_pam_sequence"].upper()
        seed = get_seed_sequence_from_protospacer(protospacer)

        if pam_rule and set(pam_rule) == {"N"}:
            full_candidates = {protospacer}
            seed_candidates = {seed}
        else:
            expanded_pams = expand_pam_iupac(pam_rule)
            full_candidates = {f"{protospacer}{pam}" for pam in expanded_pams}
            seed_candidates = {f"{seed}{pam}" for pam in expanded_pams}

        # 逆相補配列も候補に追加
        full_sequences = full_candidates | {convert_dna_to_reversed_complement_dna(seq).upper() for seq in full_candidates}
        seed_sequences = seed_candidates | {convert_dna_to_reversed_complement_dna(seq).upper() for seq in seed_candidates}

        full_sequences_per_row.append(full_sequences)
        seed_sequences_per_row.append(seed_sequences)
        all_full_sequences.update(full_sequences)
        all_seed_sequences.update(seed_sequences)

    # まず最初にAho-CorasickのAutomatonを構築
    automaton = ahocorasick.Automaton()

    for seq in all_full_sequences:
        automaton.add_word(seq, ("full", seq))

    for seq in all_seed_sequences:
        automaton.add_word(seq, ("seed", seq))

    automaton.make_automaton()

    # sgRNAごとのカウント辞書
    offtarget_count_dict_full = {seq: 0 for seq in all_full_sequences}
    offtarget_count_dict_seed = {seq: 0 for seq in all_seed_sequences}

    with open(fasta_path, 'r') as fasta_file:
        header_count = sum(1 for line in fasta_file if line.startswith(">"))
        logging.info(f"Number of chromosomes in your FASTA file: {header_count}")
        pbar = tqdm(total=header_count, desc="Calculating off-target counts", unit="chromosome")
        fasta_file.seek(0)
        chrom_seq = ""

        def process_chrom_seq(sequence):
            sequence = sequence.upper()
            # automaton.iter はマッチした箇所の (end_index, value) を返す
            for end_idx, (kind, seq) in automaton.iter(sequence):
                if kind == "full":
                    offtarget_count_dict_full[seq] += 1
                elif kind == "seed":
                    offtarget_count_dict_seed[seq] += 1

        for line in fasta_file:
            if line.startswith(">"):
                # 新しい染色体に切り替え
                if chrom_seq:
                    chrom_seq = chrom_seq.upper()
                    #automaton.iter は (end_idx, (idx, seq)) を持っていて、そこに含まれるseqが見つかったらカウントを増やす
                    process_chrom_seq(chrom_seq)
                    chrom_seq = ""
                    pbar.update(1)
            else:
                chrom_seq += line.strip()
        # 最後の染色体も処理
        if chrom_seq:
            chrom_seq = chrom_seq.upper()
            process_chrom_seq(chrom_seq)
            pbar.update(1)
            pbar.close()
    
    exploded_sgrna_df["pam+20bp_exact_match_count"] = [
        sum(offtarget_count_dict_full[sequence] for sequence in sequences)
        for sequences in full_sequences_per_row
    ]
    exploded_sgrna_df["pam+12bp_exact_match_count"] = [
        sum(offtarget_count_dict_seed[sequence] for sequence in sequences)
        for sequences in seed_sequences_per_row
    ]
    return exploded_sgrna_df

def score_offtargets(exploded_sgrna_df: pd.DataFrame, assembly_name: str, fasta_path: Path) -> pd.DataFrame:
    """
    Purpose: このモジュールのラップ関数
    """
    exploded_sgrna_df = add_crisprdirect_url_to_df(exploded_sgrna_df, assembly_name)
    exploded_sgrna_df = calculate_offtarget_site_count_ahocorasick(exploded_sgrna_df, fasta_path)
    exploded_sgrna_df = exploded_sgrna_df.drop(columns=["sgrna_target_sequence"])
    return exploded_sgrna_df
