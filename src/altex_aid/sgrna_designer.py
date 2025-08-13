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

@dataclass(frozen=True)
class BaseEditor:
    """
    BaseEditorの情報を保持するためのdataclass
    """
    base_editor_name: str # BaseEditorの名前
    pam_sequence: str # PAM配列
    editing_window_start_in_grna: int # 編集ウィンドウの開始位置 (1-indexed)
    editing_window_end_in_grna: int # 編集ウィンドウの終了位置 (1-indexed)
    base_editor_type: str # "CBE" or "ABE"

def convert_dna_to_reversed_complement_rna(sequence: str) -> str:
    """
    purpose:
        塩基配列を逆相補のRNAに変換する
    Parameters:
        sequence: 変換したいDNA配列
    Returns:
        reversed_complement_rna_sequence: 入力したDNA配列の逆相補RNA配列
    """
    complement_map = {
        "A": "U", "T": "A", "C": "G", "G": "C", "N": "N",
        "a": "u", "t": "a", "c": "g", "g": "c"
        }
    return "".join([complement_map[base] for base in reversed(sequence)])

def convert_dna_to_rna(sequence: str) -> str:
    """
    Purpose:
        塩基配列をRNAに変換する
    Parameters:
        sequence: 変換したいDNA配列
    Returns:
        rna_sequence: 入力したDNA配列のRNA配列
    """
    sequence = sequence.replace("T", "U")
    sequence = sequence.replace("t", "u")
    return sequence

def reverse_complement_pam_as_regex(pam_sequence: str) -> re.Pattern:
    """
    Purpose:
        PAM配列を逆相補鎖に変換し、Nをワイルドカードに変換
    Parameters:
        pam_sequence: str, PAM配列
    Returns:
        reversed_pam_sequence: str, 反転されたPAM配列を探すための正規表現パターン
    """
    pam_sequence = pam_sequence.upper()  # 大文字に変換
    complement_dict = {"N": "[ATGCatgc]", "A": "[Tt]", "T": "[Aa]", "G": "[Cc]", "C": "[Gg]",
                        "M": "[TtGg]", "R": "[TtCc]", "W": "[TtAa]", "S": "[CcGg]",
                        "Y": "[AaGg]", "K": "[AaCc]", "V": "[TtCcGg]", "H": "[AaTtGg]",
                        "D": "[AaTtCc]", "B": "[CcGgAa]"}
    reversed_complement_pam_regex = "".join([complement_dict[base] for base in reversed(pam_sequence)])
    return re.compile(f"(?=({reversed_complement_pam_regex}))")

def convert_pam_as_regex(pam_sequence: str) -> re.Pattern:
    """
    Purpose:
        PAM配列を正規表現パターンに変換
    Parameters:
        pam_sequence: str, PAM配列
    Returns:
        pam_regex: str, PAM配列を表す正規表現パターン
    """
    pam_sequence = pam_sequence.upper()  # 大文字に変換
    complement_dict = {"N": "[ATGCatgc]", "A": "[Aa]", "T": "[Tt]", "G": "[Gg]", "C": "[Cc]",
                        "M": "[AaCc]", "R": "[AaGg]", "W": "[AaTt]", "S": "[CcGg]",
                        "Y": "[CcTt]", "K": "[GgTt]", "V": "[AaCcGg]", "H": "[AaTtCc]",
                        "D": "[GgAaTt]", "B": "[TtGgCc]"}
    pam_regex = "".join([complement_dict[base] for base in pam_sequence])
    return re.compile(f"(?=({pam_regex}))")  


def calculate_overlap_and_unintended_edits_to_cds(
    editing_sequence: str,
    window_start_in_seq: int,
    window_end_in_seq: int,
    cds_boundary: int,
    site_type: str,
    base_editor_type: str
) -> tuple[int, int]:
    """
    Purpose:
        編集ウィンドウとCDSの重なりの長さ、およびcds中で意図しない編集を受ける可能性のある塩基の数を計算する
    Parameters:
        editing_sequence: str, 編集対象の塩基配列
        window_start_in_seq: int, 編集ウィンドウの開始位置 (0-indexed)
        window_end_in_seq: int, 編集ウィンドウの終了位置 (0-indexed)
        cds_boundary: int, CDSの境界位置 (0-indexed)
        site_type: str, "acceptor" または "donor"
    Returns:
        overlap: int, 編集ウィンドウとCDSの重なりの長さ
        unintended_edits: int, cds中の意図しない編集を受ける可能性のある塩基の数
    """
    overlap = 0
    unintended_edits = 0

    if site_type == "acceptor":
        if window_end_in_seq >= cds_boundary:
            overlap = window_end_in_seq - cds_boundary + 1
            if base_editor_type == "cbe":
                unintended_edits = editing_sequence[
                    cds_boundary: window_end_in_seq + 1  # +1はinclusiveにするため
                    ].count("G") 
            elif base_editor_type == "abe":
                unintended_edits = editing_sequence[
                    cds_boundary: window_end_in_seq + 1  # +1はinclusiveにするため
                    ].count("A")
    else:  # donor
        if window_start_in_seq <= cds_boundary:
            overlap = cds_boundary - window_start_in_seq + 1
            if base_editor_type == "cbe":
                unintended_edits = editing_sequence[
                    window_start_in_seq: cds_boundary + 1  # +1はinclusiveにするため
                ].count("G")
            elif base_editor_type == "abe":
                unintended_edits = editing_sequence[
                    window_start_in_seq: cds_boundary + 1  # +1はinclusiveにするため
                ].count("T")

    return overlap, unintended_edits

TARGET_BASE_POS = {
    ("acceptor", "cbe"): 24,
    ("donor", "cbe"): 25,
    ("acceptor", "abe"): 23,
    ("donor", "abe"): 26,
}

def decide_target_base_pos_in_sequence(
    base_editor_type: str,
    site_type:str,
) ->int:
    """
    purpose:
        base_editorの種類とサイトタイプに応じてターゲットとなる塩基の取得塩基中での位置を決定する
    Parameters:
        base_editor_type: str, "abe" または "cbe"
        site_type: str, "acceptor" または "donor"
    Returns:
        target_base_pos_in_sequence: int, 
            ターゲットとなる塩基の取得塩基中での位置 (0-indexed)
    Comments:
        取得した50塩基の配列の99%は、以下のルールに従う。
        acceptorの場合:
            23 24
        5'-- A  G  [---exon---]
        3'-- T  C  [---exon---]
        donorの場合:
                        25 26
        5'--[---exon---] G  T --3'
        3'--[---exon---] C  A --5'
    """
    return TARGET_BASE_POS[(site_type, base_editor_type.lower())]

VALID_EXON_POSITIONS = {
    "acceptor": ["internal", "last"],
    "donor": ["internal", "first"],
}

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
    return exon_position in VALID_EXON_POSITIONS[site_type]

SGRNA_START_END_RULES = {
    ("acceptor", "cbe"): lambda match: (match.end(1), match.end(1) + 20),
    ("acceptor", "abe"): lambda match: (match.start(1) - 20, match.start(1)),
    ("donor", "cbe"): lambda match: (match.end(1), match.end(1) + 20),
    ("donor", "abe"): lambda match: (match.end(1), match.end(1) + 20),
}

def decide_sgrna_start_and_end(match: re.Match, site_type: str, base_editor_type: str) -> tuple[int, int]:
    return SGRNA_START_END_RULES[(site_type, base_editor_type.lower())](match)

def design_sgrna(
    editing_sequence: str,
    pam_regex: re.Pattern,
    reversed_pam_regex: re.Pattern,
    editing_window_start_in_grna: int,
    editing_window_end_in_grna: int,
    target_base_pos_in_sequence: int,
    cds_boundary: int,
    base_editor_type: str,
    site_type: str
) -> list[SgrnaInfo]:
    """
    Purpose:
        ゲノムから取り出したSA/SD周辺配列からsgRNAの条件に適合する配列を検索する
    Parameters:
        target_sequence: str, +鎖由来の50bpのSA/SD周辺の塩基配列
        pam_sequence: str, PAM配列
        editing_window_start: int, 編集ウィンドウの開始位置(1-indexed)
        editing_window_end: int, 編集ウィンドウの終了位置(1-indexed)
        target_g_pos_in_sequence: int, 編集ターゲットとなるGの位置, 
            acceptorなら24番目のG, donorなら25番目のG (0-indexed)
        cds_boundary: int, CDSの境界位置(0-indexed) つまり、SAなら25番目の塩基、SDなら24番目の塩基
        site_type: str, "acceptor"または"donor"のどちらかを指定
    Returns:
        sgrna_list: list[SgrnaInfo], 条件に適合するsgRNAの情報のリスト
    If base_editor_type is "cbe":
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
    If base_editor_type is "abe":
        Acceptor splice siteの例:
        5'--------AG[---NGG---exon]-3'
             5'---AG----3" 20bpのsgRNA(-鎖に結合する)
        3'--------TC[---NCC---exon]-5'

        Donor splice siteの例:
        5'-[---exon--CCN---]GT--3'
                       3'---CU---5' 20bpのsgRNA(+鎖に結合する)
        3'-[---exon--GGN---]CA--5'

        この"A"をGに編集するために、+鎖に結合するsgRNAを設計する。
        つまり、+鎖を逆相補にしたものがsgRNAとなる。
        しかし、マッピング時のことを考えて、sgRNA編集ターゲット, 実際の逆相補化されたgRNA配列の両方を出力する
    """
    sgrna_list = []
    base_editor_type = base_editor_type.lower()
    if site_type == "acceptor":
        splice_site = editing_sequence[cds_boundary - 2:cds_boundary].upper()
        expected_site = "AG"
        pam_iter = reversed_pam_regex.finditer(editing_sequence) if base_editor_type == "cbe" else pam_regex.finditer(editing_sequence)
    else:
        splice_site = editing_sequence[cds_boundary + 1:cds_boundary + 3].upper()
        expected_site = "GT"
        pam_iter = reversed_pam_regex.finditer(editing_sequence)
    if splice_site != expected_site:
        return sgrna_list
    for match in pam_iter:
        sgrna_start, sgrna_end = decide_sgrna_start_and_end(match, site_type, base_editor_type)
        if sgrna_start < 0 or sgrna_end > len(editing_sequence):
            continue

        # ABEかつ、acceptorの場合だけ、sgRNAが-鎖に結合するので、+鎖上では、3'端がPAMとなる
        # そのため、編集ウィンドウの開始位置と終了位置をsgrnaの3'から数える必要がある
        if site_type == "acceptor" and base_editor_type == "abe":
            window_start_in_seq = sgrna_end - editing_window_end_in_grna
            window_end_in_seq = sgrna_end - editing_window_start_in_grna
        else:
            window_start_in_seq = sgrna_start + editing_window_start_in_grna - 1
            window_end_in_seq = sgrna_start + editing_window_end_in_grna - 1

        if not (window_start_in_seq <= target_base_pos_in_sequence <= window_end_in_seq):
            continue
        target_sequence = editing_sequence[sgrna_start:sgrna_end]
        pam_plus_target_sequence = f"{target_sequence}+{match.group(1)}" if site_type == "acceptor" and base_editor_type == "abe" else f"{match.group(1)}+{target_sequence}"
        if site_type == "acceptor" and base_editor_type == "cbe":
            actual_sequence = convert_dna_to_reversed_complement_rna(target_sequence)
        elif site_type == "acceptor" and base_editor_type == "abe":
            actual_sequence = convert_dna_to_rna(target_sequence)
        else:
            actual_sequence = convert_dna_to_reversed_complement_rna(target_sequence)
        if site_type == "acceptor" and base_editor_type == "abe":
            target_pos_in_sgrna = -(target_base_pos_in_sequence - sgrna_end)
        else:
            target_pos_in_sgrna = target_base_pos_in_sequence - sgrna_start + 1
        overlap, unintended_edits = calculate_overlap_and_unintended_edits_to_cds(
            editing_sequence=editing_sequence,
            window_start_in_seq=window_start_in_seq,
            window_end_in_seq=window_end_in_seq,
            cds_boundary=cds_boundary,
            site_type=site_type,
            base_editor_type=base_editor_type
        )
        sgrna_list.append(
            SgrnaInfo(
                target_sequence=pam_plus_target_sequence,
                actual_sequence=actual_sequence,
                start_in_sequence=sgrna_start,
                end_in_sequence=sgrna_end,
                target_pos_in_sgrna=target_pos_in_sgrna,
                overlap_between_cds_and_editing_window=overlap,
                possible_unintended_edited_base_count=unintended_edits,
            )
        )
    return sgrna_list

def design_sgrna_for_target_exon_df(
    target_exon_df: pd.DataFrame,
    pam_sequence: str,
    editing_window_start_in_grna: int,
    editing_window_end_in_grna: int,
    base_editor_type: str 
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

    pam_regex = convert_pam_as_regex(pam_sequence)
    reversed_pam_regex = reverse_complement_pam_as_regex(pam_sequence)

    def apply_design(row, site_type, base_editor_type):
        sequence_col = f"{site_type}_exon_intron_boundary_±25bp_sequence"
        if is_valid_exon_position(row["exon_position"], site_type):
            return design_sgrna(
                editing_sequence=row[sequence_col],
                pam_regex=pam_regex,
                reversed_pam_regex=reversed_pam_regex,
                editing_window_start_in_grna=editing_window_start_in_grna,
                editing_window_end_in_grna=editing_window_end_in_grna,
                target_base_pos_in_sequence=decide_target_base_pos_in_sequence(base_editor_type, site_type),
                cds_boundary=ACCEPTOR_CDS_BOUNDARY if site_type == "acceptor" else DONOR_CDS_BOUNDARY,
                base_editor_type=base_editor_type,
                site_type=site_type
            )
        return []
    
    # exontypeがa5ss-longの場合はacceptor用のsgRNAを設計しない。a5ssはacceptorの位置が-shortと同じだから。
    # exontypeがa3ss-longの場合はdonor用のsgRNAを設計しない。 a3ssはdonorの位置が-shortと同じだから。
    target_exon_df["grna_acceptor"] = target_exon_df.apply(
        lambda r:[] if r["exontype"] =="a5ss-long" else apply_design(r, "acceptor", base_editor_type), axis=1
    )
    target_exon_df["grna_donor"] = target_exon_df.apply(
        lambda r:[] if r["exontype"]=="a3ss-long" else apply_design(r, "donor", base_editor_type), axis=1
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

# 今までの動作をまとめて、sgrnaを設計する関数
# 実際に実行するのはこの関数だけでよい
def design_sgrna_for_base_editors(
    target_exon_df: pd.DataFrame,
    base_editors: list[BaseEditor],
) -> pd.DataFrame:
    """
    Purpose:
        各BaseEditorに対してsgRNAを設計し、結果をDataFrameにまとめる
    Parameters:
        target_exon_df: pd.DataFrame, 各エキソンの情報を含むDataFrame
        base_editors: list[BaseEditor], BaseEditorの情報を含むリスト
    Returns:
        pd.DataFrame, 各BaseEditorに対して設計されたsgRNAの情報を含むDataFrame
    """
    results = []  # 各BaseEditorの結果を格納するリスト

    for base_editor in base_editors:
        # 1. sgRNAを設計する
        temp_df = design_sgrna_for_target_exon_df(
            target_exon_df=target_exon_df,
            pam_sequence=base_editor.pam_sequence,
            editing_window_start_in_grna=base_editor.editing_window_start_in_grna,
            editing_window_end_in_grna=base_editor.editing_window_end_in_grna,
            base_editor_type=base_editor.base_editor_type
        )
        # 2. sgRNAの情報を展開する
        temp_df = organize_target_exon_df_with_grna_sequence(temp_df)
        # 3. sgRNAの開始位置と終了位置をゲノム上の位置に変換する
        temp_df = convert_sgrna_start_end_position_to_position_in_chromosome(temp_df)

        # 4. 列名にbase_editorの名前を付ける
        temp_df = temp_df.rename(
            columns={
                col: f"{base_editor.base_editor_name}_{col}" 
                for col in temp_df.columns
                if col.startswith("acceptor_sgrna") or col.startswith("donor_sgrna")
            }
        )

        # 結果をリストに追加
        results.append(temp_df)

    # すべての結果を結合
    final_result = pd.concat(results, axis=1)
    final_result = final_result.loc[:, ~final_result.columns.duplicated()]
    return final_result.reset_index(drop=True)

# 論文用には、各BEのsgRNAを1行にまとめたほうが便利ではあるが、ユーザーにはそれは必要ない。のちのexplodeを考えて、main.pyではdictで返すようにする
def design_sgrna_for_base_editors_dict(
    target_exon_df: pd.DataFrame,
    base_editors: list[BaseEditor],
) -> dict[str, pd.DataFrame]:
    """
    Purpose:
        各BaseEditorに対してsgRNAを設計し、結果をdictにまとめる
    Parameters:
        target_exon_df: pd.DataFrame, 各エキソンの情報を含むDataFrame
        base_editors: list[BaseEditor], BaseEditorの情報を含むリスト
    Returns:
        dict[str, pd.DataFrame], 各BaseEditorに対して設計されたsgRNAの情報を含むDataFrame
    """
    results = {}  # 各BaseEditorの結果を格納する辞書

    for base_editor in base_editors:
        # 1. sgRNAを設計する
        temp_df = design_sgrna_for_target_exon_df(
            target_exon_df=target_exon_df,
            pam_sequence=base_editor.pam_sequence,
            editing_window_start_in_grna=base_editor.editing_window_start_in_grna,
            editing_window_end_in_grna=base_editor.editing_window_end_in_grna,
            base_editor_type=base_editor.base_editor_type
        )
        # 2. sgRNAの情報を展開する
        temp_df = organize_target_exon_df_with_grna_sequence(temp_df)
        # 3. sgRNAの開始位置と終了位置をゲノム上の位置に変換する
        temp_df = convert_sgrna_start_end_position_to_position_in_chromosome(temp_df)

        # 結果を辞書に追加
        results[base_editor.base_editor_name] = temp_df

    return results


def make_preset_base_editors() -> dict[str, BaseEditor]:
    """
    Purpose:
        デフォルトのBaseEditorのリストを返す
    Returns:
        dict[str, BaseEditor], デフォルトのBaseEditorのリスト
    """
    return {
        "target_aid": BaseEditor(
            base_editor_name="target_aid",
            pam_sequence="NGG",
            editing_window_start_in_grna=17,
            editing_window_end_in_grna=19,
            base_editor_type="cbe"
        ),
        "be4max": BaseEditor(
            base_editor_name="be4max",
            pam_sequence="NGG",
            editing_window_start_in_grna=12,
            editing_window_end_in_grna=17,
            base_editor_type="cbe"
        ),
        "abe8e": BaseEditor(
            base_editor_name="abe8e",
            pam_sequence="NGG",
            editing_window_start_in_grna=12,
            editing_window_end_in_grna=17,
            base_editor_type="abe"
        ),
    }