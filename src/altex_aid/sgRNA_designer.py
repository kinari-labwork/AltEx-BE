from __future__ import annotations
import pandas as pd
import re

def get_reversed_complement(sequence: str) -> str:
    """
    Purpose:
        塩基配列の逆相補鎖を取得する
    Parameters:
        sequence: str, 塩基配列
    Returns:
        reversed_complement: str, 逆相補鎖の塩基配列
    """
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    reversed_complement = "".join([complement_map[base] for base in reversed(sequence)])
    return reversed_complement

def reverse_pam_sequence(pam_sequence:str)->str:
    """
    Purpose:
        PAM配列を逆相補鎖に変換し、Nをワイルドカードに変換
    Parameters:
        pam_sequence: str, PAM配列
    Returns:
        reversed_pam_sequence: str, 反転されたPAM配列
    """
    complement_dict = {"N":"[ATGC]","A":"[Tt]", "T":"[Aa]", "G":"[Cc]", "C":"[Gg]"}
    reversed_pam_sequence = "".join([complement_dict[base] for base in reversed(pam_sequence)])
    return reversed_pam_sequence

def design_sgrna_for_acceptor_with_cbe(
    editing_sequence: str,
    pam_sequence: str,
    editing_window_start_in_grna: int,
    editing_window_end_in_grna: int,
) -> list[tuple[str, str, int, int, int]]:
    """
    Purpose:
        ゲノムから取り出したSA周辺配列からsgRNAの条件に適合する配列を検索する
    Parameters:
        target_sequence: str, +鎖由来の50bpのSA/SD周辺の塩基配列
        pam_sequence: str, PAM配列
        editing_window_start: int, 編集ウィンドウの開始位置(0-indexed)
        editing_window_end: int, 編集ウィンドウの終了位置(0-indexed)
    Returns:
        sgRNAs: tupleのlist, 各要素は(sgrnaのターゲットとなる塩基配列, sgrnaの実際の配列（逆相補にしたもの）, 開始位置, 終了位置, cdsと重なる塩基の数)の形式
        開始・終了位置は取得した配列内での位置なので、後で正確な位置に変換する必要がある
    Comments:
        5'---CCN---AG[---exon]-3'
                -------- 20bpのsgRNA(+鎖に結合する)
        3'---GGN---TC[---exon]-5'
        この"C"をTに編集するために、+鎖に結合するsgRNAを設計するべき。
        つまり、+鎖を逆相補にしたものがsgRNAとなる。
        しかし、マッピング時のことを考えて、sgRNA編集ターゲット, 実際の逆相補化されたgRNA配列の両方を出力する
    """

    # 変換対象は＋鎖のGT-AGルールには存在しないので、-鎖のCをTに編集する必要がある
    # 編集対象の塩基が存在する-鎖にPAM配列が存在する必要があるため、+鎖ではPAMの逆相補的な塩基を探す必要がある
    reversed_pam_sequence = reverse_pam_sequence(pam_sequence)

    sgrna_list = [] #候補が全く見つからなかった場合は空のリストを返す
    target_pos_g = 24 # 編集対象の塩基位置（0-indexed）
    cds_start_index = 25 # CDSの開始位置（0-indexed）

    if editing_sequence[23:25].upper() == "AG":
        for match in re.finditer(reversed_pam_sequence, editing_sequence):
            grna_candidate_start_in_sequence = match.end()
            grna_candidate_end_in_sequence = grna_candidate_start_in_sequence + 20

            # grna_candidate_start_indexがediting_sequenceの長さを超える場合はスキップ
            if grna_candidate_end_in_sequence > len(editing_sequence):
                continue

            # 編集対象のGの位置がgRNA候補内にない場合はスキップ
            if not (grna_candidate_start_in_sequence <= target_pos_g < grna_candidate_end_in_sequence):
                continue

            # 編集ウィンドウの開始位置と終了位置をediting_sequence内での位置に変換
            editing_window_start_in_sequence = grna_candidate_start_in_sequence + editing_window_start_in_grna - 1
            editing_window_end_in_sequence = grna_candidate_start_in_sequence + editing_window_end_in_grna - 1

            # 編集ウィンドウ内に編集対象のGがあるか確認
            if editing_window_start_in_sequence <= target_pos_g <= editing_window_end_in_sequence:
                sgrna_target_sequence = editing_sequence[grna_candidate_start_in_sequence : grna_candidate_end_in_sequence]
                
                #sgRNAは取得したプラス鎖の逆相補鎖となる
                actual_sgrna_sequence = get_reversed_complement(sgrna_target_sequence)
                
                #編集ウインドウとCDSの重なりを計算
                if editing_window_start_in_sequence > cds_start_index :
                    overlap_between_cds_and_editing_window = 0
                else:
                    overlap_between_cds_and_editing_window = cds_start_index - editing_window_start_in_sequence +1
                    #重なりがある場合、そこに含まれるGの数をカウントする
                    possible_unintended_edited_base_count = editing_sequence[cds_start_index:cds_start_index + overlap_between_cds_and_editing_window].count('G')

                sgrna_list.append((
                    sgrna_target_sequence,
                    actual_sgrna_sequence, 
                    grna_candidate_start_in_sequence, 
                    grna_candidate_end_in_sequence, 
                    overlap_between_cds_and_editing_window,
                    possible_unintended_edited_base_count
                    ))

    return sgrna_list




def design_sgrna_for_donor_with_cbe(
    editing_sequence: str,
    pam_sequence: str,
    editing_window_start_in_grna: int,
    editing_window_end_in_grna: int,
    ) -> list[tuple[str, str, int, int, int, int, int]]:
    """
    Purpose:
        ゲノムから取り出したSD周辺配列からsgRNAの条件に適合する配列を検索する
    Parameters:
        target_sequence: str, +鎖由来の50bpのSA/SD周辺の塩基配列
        pam_sequence: str, PAM配列
        editing_window_start: int, 編集ウィンドウの開始位置(0-indexed)
        editing_window_end: int, 編集ウィンドウの終了位置(0-indexed)
    Returns:
        sgRNAs: tupleのlist, 各要素は(sgrnaのターゲットとなる塩基配列, 実際のsgrna配列（逆相補にしたもの）,開始位置, 終了位置,cdsと重なる塩基の数, 意図しないcds内での変異が起こりうる塩基の数)の形式
        開始・終了位置は取得した配列内での位置なので、後で正確な位置に変換する必要がある
    Comments:
        5'-[---exon--CCN---]GT--3'
                        -------- 20bpのsgRNA(+鎖に結合する)
        3'-[---exon--GGN---]CA--5'
        この"C"をTに編集するために、+鎖に結合するsgRNAを設計する,
        つまり、+鎖を逆相補にしたものがsgRNAとなる
    """

    # 変換対象は＋鎖のGT-AGルールには存在しないので、-鎖のCをTに編集する必要がある
    # 編集対象の塩基が存在する-鎖にPAM配列が存在する必要があるため、+鎖ではPAMの逆相補的な塩基を探す必要がある
    reversed_pam_sequence = reverse_pam_sequence(pam_sequence)

    sgRNA = []
    target_pos_g = 25 # 編集対象の塩基位置（0-indexed）
    cds_end_index = 24 # CDSの開始位置（0-indexed）

    if editing_sequence[25:27].upper() == "GT":
        for match in re.finditer(reversed_pam_sequence, editing_sequence):
            grna_candidate_start_in_sequence = match.end()
            grna_candidate_end_in_sequence = grna_candidate_start_in_sequence + 20

            # grna_candidate_start_indexがediting_sequenceの長さを超える場合はスキップ
            if grna_candidate_end_in_sequence > len(editing_sequence):
                continue

            # 編集対象のGの位置がgRNA候補内にない場合はスキップ
            if not (grna_candidate_start_in_sequence <= target_pos_g < grna_candidate_end_in_sequence):
                continue

            # 編集ウィンドウの開始位置と終了位置をediting_sequence内での位置に変換
            editing_window_start_in_sequence = grna_candidate_start_in_sequence + editing_window_start_in_grna - 1
            editing_window_end_in_sequence = grna_candidate_start_in_sequence + editing_window_end_in_grna - 1

            # 編集ウィンドウ内に編集対象のGがあるか確認
            if editing_window_start_in_sequence <= target_pos_g <= editing_window_end_in_sequence:
                sgRNA_sequence = editing_sequence[grna_candidate_start_in_sequence : grna_candidate_end_in_sequence]
                
                #sgRNAは取得したプラス鎖の逆相補鎖となるため、変換する
                sgRNA_sequence = get_reversed_complement(sgRNA_sequence)
                
                #sgrnaの何塩基目に編集対象塩基が来るかを計算
                target_pos_g_in_sgrna = target_pos_g - grna_candidate_start_in_sequence + 1

                # sgRNAとCDSの重なりを計算
                if editing_window_start_in_sequence > cds_end_index :
                    overlap_between_cds_and_editing_window = 0
                else:
                    overlap_between_cds_and_editing_window = cds_end_index - editing_window_start_in_sequence +1

                if editing_window_start_in_sequence > cds_end_index :
                    overlap_between_cds_and_editing_window = 0
                else:
                    overlap_between_cds_and_editing_window = cds_end_index - editing_window_start_in_sequence +1
                    #重なりがある場合、そこに含まれるGの数をカウントする
                    possible_unintended_edited_base_count = editing_sequence[cds_end_index - overlap_between_cds_and_editing_window +1 :cds_end_index +1].count('G')                
                
                sgRNA.append((
                    sgRNA_sequence, 
                    grna_candidate_start_in_sequence, 
                    grna_candidate_end_in_sequence, 
                    target_pos_g_in_sgrna,
                    overlap_between_cds_and_editing_window,
                    possible_unintended_edited_base_count
                    ))

    return sgRNA

def design_sgrna_for_target_exon_df_with_sequence(
    target_exon_df_with_acceptor_and_donor_sequence: pd.DataFrame,
    pam_sequence: str,
    editing_window_start_in_grna: int,
    editing_window_end_in_grna: int,
) -> pd.DataFrame:
    """
    Purpose:
        target_exon_dfの各行に対してsgRNAを設計する
    Parameters:
        target_exon_df: pd.DataFrame, ターゲットエクソンの情報を含むDataFrame
        pam_sequence: str, PAM配列
        editing_window_start_in_grna: int, 編集ウィンドウの開始位置(0-indexed)
        editing_window_end_in_grna: int, 編集ウィンドウの終了位置(0-indexed)
    Returns:
        pd.DataFrame, sgRNAの情報を含むDataFrame
    """
    
    target_exon_df_with_acceptor_and_donor_sequence["grna_acceptor"] = target_exon_df_with_acceptor_and_donor_sequence.apply(
        lambda row: design_sgrna_for_acceptor_with_cbe(
            editing_sequence=row["acceptor_sequence"],
            pam_sequence=pam_sequence,
            editing_window_start_in_grna=editing_window_start_in_grna,
            editing_window_end_in_grna=editing_window_end_in_grna)
            if row['exon_position'] in ["internal", "last", "a3ss-long"]
            else [],
        axis=1,
        )
    
    target_exon_df_with_acceptor_and_donor_sequence["grna_donor"] = target_exon_df_with_acceptor_and_donor_sequence.apply(
        lambda row: design_sgrna_for_donor_with_cbe(
            editing_sequence=row["donor_sequence"],
            pam_sequence=pam_sequence,
            editing_window_start_in_grna=editing_window_start_in_grna,
            editing_window_end_in_grna=editing_window_end_in_grna)
            if row['exon_position'] in ["internal", "first", "a5ss-long"]
            else [],
        axis=1,
        )

    target_exon_df_with_grna_sequence = target_exon_df_with_acceptor_and_donor_sequence

    return target_exon_df_with_grna_sequence.reset_index(drop=True)

def extract_acceptor_sgrna_features(sgrna_list):
    """acceptor用sgRNAリストから各特徴を抽出"""
    if not sgrna_list:
        return [], [], [], [], [], []
    return (
        [t[0] for t in sgrna_list],  # target_sequence
        [t[1] for t in sgrna_list],  # actual_sequence
        [t[2] for t in sgrna_list],  # start_in_sequence
        [t[3] for t in sgrna_list],  # end_in_sequence
        [t[4] for t in sgrna_list],  # overlap_between_cds_and_editing_window
        [t[5] for t in sgrna_list],  # possible_unintended_edited_base_count
    )

def extract_donor_sgrna_features(sgrna_list):
    """donor用sgRNAリストから各特徴を抽出"""
    if not sgrna_list:
        return [], [], [], [], [], []
    return (
        [t[0] for t in sgrna_list],  # target_sequence
        [t[1] for t in sgrna_list],  # start_in_sequence
        [t[2] for t in sgrna_list],  # end_in_sequence
        [t[3] for t in sgrna_list],  # target_pos_in_sgrna
        [t[4] for t in sgrna_list],  # overlap_between_cds_and_editing_window
        [t[5] for t in sgrna_list],  # possible_unintended_edited_base_count
    )


def organize_target_exon_df_with_grna_sequence(
    target_exon_df_with_grna_sequence: pd.DataFrame,
) -> pd.DataFrame:
    """grna列から各特徴量列を展開"""
    (
        target_exon_df_with_grna_sequence["acceptor_sgrna_target_sequence"],
        target_exon_df_with_grna_sequence["acceptor_sgrna_actual_sequence"],
        target_exon_df_with_grna_sequence["acceptor_sgrna_start_in_sequence"],
        target_exon_df_with_grna_sequence["acceptor_sgrna_end_in_sequence"],
        target_exon_df_with_grna_sequence["acceptor_sgrna_overlap_between_cds_and_editing_window"],
        target_exon_df_with_grna_sequence["acceptor_sgrna_possible_unintended_edited_base_count"],
    ) = zip(*target_exon_df_with_grna_sequence["grna_acceptor"].apply(extract_acceptor_sgrna_features))
    (
        target_exon_df_with_grna_sequence["donor_sgrna_target_sequence"],
        target_exon_df_with_grna_sequence["donor_sgrna_start_in_sequence"],
        target_exon_df_with_grna_sequence["donor_sgrna_end_in_sequence"],
        target_exon_df_with_grna_sequence["donor_sgrna_target_pos_in_sgrna"],
        target_exon_df_with_grna_sequence["donor_sgrna_overlap_between_cds_and_editing_window"],
        target_exon_df_with_grna_sequence["donor_sgrna_possible_unintended_edited_base_count"],
    ) = zip(*target_exon_df_with_grna_sequence["grna_donor"].apply(extract_donor_sgrna_features))
    # 不要な列を削除
    return target_exon_df_with_grna_sequence.drop(columns=["grna_acceptor", "grna_donor"]).reset_index(drop=True)


def modify_sgrna_start_end_position(
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
    target_exon_df_with_grna_sequence["acceptor_sgrna_start_in_genome"] = (
        target_exon_df_with_grna_sequence["acceptor_sgrna_start_in_sequence"]
        + target_exon_df_with_grna_sequence["chromStart_acceptor"]
    )
    target_exon_df_with_grna_sequence["acceptor_sgrna_end_in_genome"] = (
        target_exon_df_with_grna_sequence["chromStart_acceptor"]
        + target_exon_df_with_grna_sequence["acceptor_sgrna_end_in_sequence"]
        
    )
    
    target_exon_df_with_grna_sequence["donor_sgrna_start_in_genome"] = (
        target_exon_df_with_grna_sequence["donor_sgrna_start_in_sequence"]
        + target_exon_df_with_grna_sequence["chromStart_donor"]

    )
    target_exon_df_with_grna_sequence["donor_sgrna_end_in_genome"] = (
        target_exon_df_with_grna_sequence["chromEnd_donor"]
        + target_exon_df_with_grna_sequence["donor_sgrna_end_in_sequence"]
    )

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

