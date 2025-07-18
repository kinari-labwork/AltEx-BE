from __future__ import annotations
from dataclasses import (dataclass,astuple)
import pandas as pd
import re


@dataclass(frozen=True)
class SgrnaInfo:
    """
    使うデータをdataclassで定義することで、コードの可読性を向上させる。
    """

    target_sequence: str # sgRNAのターゲットとなる塩基配列
    actual_sequence: str # 実際に設計されるsgRNAの配列 (+鎖に対して結合するため、逆相補化される)
    start_in_sequence: int # 取得しているSA/SD周辺の50塩基対のうち、sgRNAの開始位置 (0-indexed)
    end_in_sequence: int # 取得しているSA/SD周辺の50塩基対のうち、sgRNAの終了位置 (0-indexed)
    target_pos_in_sgrna: int # sgRNAの中で、編集ターゲットとなる塩基の位置 (1-indexed)
    overlap_between_cds_and_editing_window: int # 編集ウィンドウとCDSの重なりの長さ
    possible_unintended_edited_base_count: int # 意図しないcdsでの変異が起こる可能性のある塩基の数


def get_reversed_complement(sequence: str) -> str:
    """
    purpose:
        塩基配列を逆相補のRNAに変換する
    Parameters:
        sequence: 変換したい塩基配列
    Returns:
        逆相補の塩基配列
    """
    complement_map = {
        "A": "U", "T": "A", "C": "G", "G": "C", "N": "N",
        "a": "u", "t": "a", "c": "g", "g": "c"
        }
    return "".join([complement_map[base] for base in reversed(sequence)])


def reverse_complement_pam_as_regex(pam_sequence: str) -> str:
    """
    Purpose:
        PAM配列を逆相補鎖に変換し、Nをワイルドカードに変換
    Parameters:
        pam_sequence: str, PAM配列
    Returns:
        reversed_pam_sequence: str, 反転されたPAM配列を探すための正規表現パターン
    """
    complement_dict = {"N": "[ATGCatgc]", "A": "[Tt]", "T": "[Aa]", "G": "[Cc]", "C": "[Gg]"}
    return "".join([complement_dict[base] for base in reversed(pam_sequence)])


def design_sgrna(
    editing_sequence: str,
    pam_sequence: str,
    editing_window_start_in_grna: int, # 1-indexed
    editing_window_end_in_grna: int, # 1-indexed
    target_g_pos_in_sequence: int,
    cds_boundary: int,
    site_type: str,
) -> list[SgrnaInfo]:
    """
    Purpose:
        ゲノムから取り出したSA/SD周辺配列からsgRNAの条件に適合する配列を検索する
    Parameters:
        target_sequence: str, +鎖由来の50bpのSA/SD周辺の塩基配列
        pam_sequence: str, PAM配列
        editing_window_start: int, 編集ウィンドウの開始位置(1-indexed)
        editing_window_end: int, 編集ウィンドウの終了位置(1-indexed)
        splice_site_pos: int, 5'から数えたときのSA/SDの開始位置(0-indexed) つまり、SAなら23番目のA, SDなら25番目のG
        cds_boundary: int, CDSの境界位置(0-indexed) つまり、SAなら25番目の塩基、SDなら24番目の塩基
        site_type: str, "acceptor"または"donor"のどちらかを指定
    Returns:
        sgrna_list: list[SgrnaInfo], 条件に適合するsgRNAの情報のリスト
    Comments:
        Acceptor splice siteの例:
        5'---CCN---AG[---exon]-3'
              3'---UC---5" 20bpのsgRNA(+鎖に結合する)
        3'---GGN---TC[---exon]-5'

        Donor splice siteの例:
        5'-[---exon--CCN---]GT--3'
                       3'---CU---5' 20bpのsgRNA(+鎖に結合する)
        3'-[---exon--GGN---]CA--5'

        この"C"をTに編集するために、+鎖に結合するsgRNAを設計する。
        つまり、+鎖を逆相補にしたものがsgRNAとなる。
        しかし、マッピング時のことを考えて、sgRNA編集ターゲット, 実際の逆相補化されたgRNA配列の両方を出力する
    """
    reversed_pam = reverse_complement_pam_as_regex(pam_sequence)
    reversed_pam = f"(?=({reversed_pam}))"  # lookaheadを使って、重複を許して探索する
    sgrna_list = []

    splice_site = editing_sequence[target_g_pos_in_sequence - 1: target_g_pos_in_sequence + 1].upper() if site_type == "acceptor" else editing_sequence[target_g_pos_in_sequence : target_g_pos_in_sequence + 2].upper()
    expected_site = "AG" if site_type == "acceptor" else "GT"
    if splice_site != expected_site:
        return sgrna_list

    for match in re.finditer(reversed_pam, editing_sequence):
        grna_start = match.end(1) 
        grna_end = grna_start + 20
        if grna_end > len(editing_sequence):
            continue

        if not (grna_start <= target_g_pos_in_sequence < grna_end):
            continue

        window_start_in_seq = grna_start + editing_window_start_in_grna -1 # 1-indexedから0-indexedに変換するため、余分に1を引く
        window_end_in_seq = grna_start + editing_window_end_in_grna -1
        if not (window_start_in_seq <= target_g_pos_in_sequence <= window_end_in_seq):
            continue

        target_sequence = editing_sequence[grna_start:grna_end]
        actual_sequence = get_reversed_complement(target_sequence)
        target_pos_in_sgrna = target_g_pos_in_sequence - grna_start + 1

        overlap = 0
        unintended_edits = 0
        if site_type == "acceptor":
            if window_end_in_seq >= cds_boundary:
                overlap = window_end_in_seq - cds_boundary +1
                unintended_edits = editing_sequence[
                    cds_boundary: window_end_in_seq +1 # +1はinclusiveにするため
                ].count("G")
        else:  # donor
            if window_start_in_seq <= cds_boundary:
                overlap = cds_boundary - window_start_in_seq +1
                unintended_edits = editing_sequence[
                    window_start_in_seq: cds_boundary + 1 # +1はinclusiveにするため
                ].count("G")

        sgrna_list.append(
            SgrnaInfo(
                target_sequence=target_sequence,
                actual_sequence=actual_sequence,
                start_in_sequence=grna_start,
                end_in_sequence=grna_end,
                target_pos_in_sgrna=target_pos_in_sgrna,
                overlap_between_cds_and_editing_window=overlap,
                possible_unintended_edited_base_count=unintended_edits,
            )
        )
        
    return sgrna_list


def is_valid_exon_position(exon_position: str, site_type: str) -> bool:
    """
    Purpose:
        exon_position が site_type に応じて有効かどうかを判定する
    Parameters:
        exon_position: str, エキソンの位置 ("internal", "first", "last" のいずれか)
        site_type: str, "acceptor" または "donor"
    Returns:
        bool, 有効なら True, 無効なら False
    """
    valid_positions = ["internal", "last"] if site_type == "acceptor" else ["internal", "first"]
    return exon_position in valid_positions



def design_sgrna_for_target_exon_df(
    target_exon_df: pd.DataFrame,
    pam_sequence: str,
    editing_window_start_in_grna: int,
    editing_window_end_in_grna: int,
) -> pd.DataFrame:
    """
    Purpose:
        各エキソンのSA/SD周辺配列が格納されているDataFrameに対して、sgRNAを設計する
    Parameters:
        target_exon_df: DataFrame, 各エキソンの情報を含むDataFrame
        pam_sequence: str, PAM配列
        editing_window_start_in_grna:編集ウィンドウの開始位置(1-indexed)
        editing_window_end_in_grna: 編集ウィンドウの終了位置(1-indexed)
    Returns:
        各エキソンに対して設計されたsgRNAの情報を含むDataFrame
    """
    ACCEPTOR_CDS_BOUNDARY = 25 # 25番目以後の塩基がCDSに含まれる
    DONOR_CDS_BOUNDARY = 24 # 24番目以前の塩基がCDSに含まれる

    def apply_design(row, site_type):
        sequence_col = f"{site_type}_sequence"
        if is_valid_exon_position(row["exon_position"], site_type):
            return design_sgrna(
                editing_sequence=row[sequence_col],
                pam_sequence=pam_sequence,
                editing_window_start_in_grna=editing_window_start_in_grna,
                editing_window_end_in_grna=editing_window_end_in_grna,
                target_g_pos_in_sequence=(
                    24 if site_type == "acceptor" else 25 # acceptorなら24番目のG, donorなら25番目のGが編集ターゲット
                ),
                cds_boundary=(
                    ACCEPTOR_CDS_BOUNDARY
                    if site_type == "acceptor"
                    else DONOR_CDS_BOUNDARY
                ),
                site_type=site_type,
            )
        return []
    
    # exontypeがa5ss-longの場合はacceptor用のsgRNAを設計しない。a5ssはacceptorの位置が-shortと同じだから。
    # exontypeがa3ss-longの場合はdonor用のsgRNAを設計しない。 a3ssはdonorの位置が-shortと同じだから。
    target_exon_df["grna_acceptor"] = target_exon_df.apply(
        lambda r:[] if r["exontype"] =="a5ss-long" else apply_design(r, "acceptor"), axis=1
    )
    target_exon_df["grna_donor"] = target_exon_df.apply(
        lambda r:[] if r["exontype"]=="a3ss-long" else apply_design(r, "donor"), axis=1
    )
    return target_exon_df




def extract_sgrna_features(sgrna_list: list[SgrnaInfo]) -> tuple[list,list,list,list,list,list,list]:
    """
    Purpose:
    donor用sgRNAリストから各特徴を抽出
    Parameters:
        sgrna_list: list[SgrnaInfo], donor用sgRNAinfoのリスト
    Returns:
        各sgRNAの特徴を抽出したリスト
        各リストの要素は、sgRNAのターゲット配列、
        実際の配列、開始位置、終了位置、ターゲット位置、
        編集ウィンドウとCDSの重なり、意図しない編集塩基の数
    """
    if not sgrna_list:
        return tuple([] for _ in range(7))  # 空のリストを返す
    return tuple(map(list,zip(*map(astuple, sgrna_list))))

def organize_target_exon_df_with_grna_sequence(target_exon_df_with_grna_sequence: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
    grna列から各特徴量列を展開する（1エキソン1行、各列はリスト）
    Parameters:
        target_exon_df_with_grna_sequence: pd.DataFrame, sgRNAの情報を含むDataFrame
    Returns:
        pd.DataFrame, 各エキソンに対してsgRNAの情報を展開したDataFrame
    Comments:
        各エキソンに対して、acceptorとdonorのsgRNA情報を展開する。
        各列はリストで、sgRNAのターゲット配列、実際の配列、開始位置、終了位置、ターゲット位置
        編集ウィンドウとCDSの重なり、意図しない編集塩基の数を含む。
    """
    for site in ["acceptor", "donor"]:
        extracted = zip(*target_exon_df_with_grna_sequence[f"grna_{site}"].apply(extract_sgrna_features))
        columns =[
            f"{site}_sgrna_target_sequence",
            f"{site}_sgrna_actual_sequence",
            f"{site}_sgrna_start_in_sequence",
            f"{site}_sgrna_end_in_sequence",
            f"{site}_sgrna_target_pos_in_sgrna",
            f"{site}_sgrna_overlap_between_cds_and_editing_window",
            f"{site}_sgrna_possible_unintended_edited_base_count",
        ]
        for col, values in zip(columns, extracted):
            target_exon_df_with_grna_sequence[col] = values
    # 不要な列を削除
    return target_exon_df_with_grna_sequence.drop(columns=["grna_acceptor", "grna_donor"]).reset_index(drop=True)

def convert_sgrna_start_end_position_to_position_in_chromosome(
    target_exon_df_with_grna_sequence: pd.DataFrame,
) -> pd.DataFrame:
    """
    Purpose:
        この段階でのacceptor/donor_start_in_sequenceは、ゲノム上の位置ではなく、取得した配列内での位置である。
        sgRNAの開始位置と終了位置を、ゲノム上の位置に変換する
    Parameters:
        target_exon_df_with_grna_sequence: pd.DataFrame, sgRNAの情報を含むDataFrame
    Returns:
        pd.DataFrame, sgRNAの開始位置と終了位置がゲノム上の位置に変換されたDataFrame
    """
    for splicesite in ["acceptor", "donor"]:
        target_exon_df_with_grna_sequence[f"{splicesite}_sgrna_start_in_genome"] = [
            [start + chrom_start for start in starts]
        for starts, chrom_start in zip(
            target_exon_df_with_grna_sequence[f"{splicesite}_sgrna_start_in_sequence"],
            target_exon_df_with_grna_sequence[f"chromStart_{splicesite}"]
        )
    ]
        target_exon_df_with_grna_sequence[f"{splicesite}_sgrna_end_in_genome"] = [
            [end + chrom_start for end in ends]
        for ends, chrom_start in zip(
            target_exon_df_with_grna_sequence[f"{splicesite}_sgrna_end_in_sequence"],
            target_exon_df_with_grna_sequence[f"chromStart_{splicesite}"]
        )
    ]

    # 不要な列を削除
    target_exon_df_with_grna_sequence = target_exon_df_with_grna_sequence.drop(
        columns=[
            "acceptor_sgrna_start_in_sequence",
            "acceptor_sgrna_end_in_sequence",
            "donor_sgrna_start_in_sequence",
            "donor_sgrna_end_in_sequence",
            "chromStart_acceptor",
            "chromEnd_acceptor",
            "chromStart_donor",
            "chromEnd_donor",
        ]
    )

    return target_exon_df_with_grna_sequence.reset_index(drop=True)