import pandas as pd

def prioritize_sgrna(exploded_sgrna_with_offtarget_info: pd.DataFrame) -> pd.DataFrame:
    """
    Exon ごとに sgRNA を優先順位付けする。
    
    ソート順序（優先順位が高い順）:
    1. pam+20bp_exact_match_count（特異的な完全一致、小さいほど優先）
    2. gc_valid（GC% が 40-60% 範囲内か、範囲内が優先）
    3. pam+12bp_exact_match_count（より緩いマッチ、小さいほど優先）
    4. sgrna_possible_unintended_edited_base_count（CDS内の予期しない編集、タイの場合のみ）
    
    Parameters:
        exploded_sgrna_with_offtarget_info: sgRNA データフレーム
        
    Returns:
        優先度を付けたデータフレーム
    """
    df = exploded_sgrna_with_offtarget_info.copy()
    
    # GC content が最適範囲（40-60%）かどうか（バイナリ判定）
    gc_content = df["sgrna_sequence"].apply(
        lambda seq: (seq.upper().count('G') + seq.upper().count('C')) / len(seq) * 100
        if len(seq) > 0 else 0
    )
    df["gc_valid"] = (gc_content >= 40) & (gc_content <= 60)
    
    # 5. Exon ごとにソート
    sorted_sgrna_df = df.sort_values(
        by=[
            "geneName",
            "exonStarts",
            "exonEnds",
            "pam+20bp_exact_match_count",  # オフターゲットの20bp完全一致数（小さいほど良い）
            "gc_valid",  # GC content が 40-60% 範囲内か（True が優先）
            "pam+12bp_exact_match_count", # オフターゲットの12bp完全一致数（小さいほど良い）
            "sgrna_possible_unintended_edited_base_count",  # タイの時の補助基準
        ],
        ascending=[True, True, True, True, False, True, True],  # gc_valid は False > True なので False に
    ).reset_index(drop=True)
    
    # 6. Exon ごとに優先度を付与
    sorted_sgrna_df["sgrna_priority"] = (
        sorted_sgrna_df
        .groupby(["geneName", "exonStarts", "exonEnds"])
        .cumcount() + 1
    )
    
    # 7. 中間計算カラムを削除
    sorted_sgrna_df = sorted_sgrna_df.drop(
        columns=["gc_valid"]
    )
    
    return sorted_sgrna_df
