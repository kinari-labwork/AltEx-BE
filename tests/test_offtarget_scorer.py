import pandas as pd
import numpy as np
from altex_aid.offtarget_scorer import (
    add_crisprdirect_url_to_df,
    calculate_offtarget_site_count
)

def test_add_crisprdirect_url_to_df_robust_success():
    # 1. Arrange: テストデータの準備
    input_df = pd.DataFrame({
        "sgrna_target_sequence": ["GAT+TACA", "ATA+TATA"],
        "base_editor_pam": ["NGG", "NGG"],
    })
    
    expected_df = pd.DataFrame({
        "sgrna_target_sequence": ["GAT+TACA", "ATA+TATA"],
        "base_editor_pam": ["NGG", "NGG"],
        "crisprdirect_url": [
            "https://crispr.dbcls.jp/?userseq=gattaca&pam=NGG&db=hg38",
            "https://crispr.dbcls.jp/?userseq=atatata&pam=NGG&db=hg38",
        ]
    })

    # 2. Act: テスト対象の関数を実行
    output_df = add_crisprdirect_url_to_df(input_df, assembly_name="hg38")

    # 3. Assert: 結果を検証
    pd.testing.assert_frame_equal(output_df, expected_df)


def test_calculate_offtarget_site_count():
    """
    ダミーのtest.fastaの中身
    >chr1
    AGCTAGCTAGCTAGCTAGCT
    >chr2
    GATTACAGATTACAGATTACA
    """
    input_df = pd.DataFrame({
        "uuid": ["id_01", "id_02", "id_03", "id_04", "id_05"],
        "sgrna_target_sequence": [
            "GAT+TACA",      # test.fa に3回出現
            "A+GCT",         # test.fa に5回出現
            "TTT+TTTT",      # test.fa には出現しない (0回)
            "GAT+TACA",      # 重複ケース
            pd.NA           # 欠損値ケース
        ]
    })
    expected_df = pd.DataFrame({
    "uuid": ["id_01", "id_02", "id_03", "id_04", "id_05"],
    "sgrna_target_sequence": [
        "GAT+TACA", "A+GCT", "TTT+TTTT", "GAT+TACA", pd.NA
    ],
    "pam+20bp_exact_match_count": [3.0, 5.0, 0.0, 3.0, np.nan]
    })
    output_df = calculate_offtarget_site_count(input_df, "tests/data/test.fasta")
    pd.testing.assert_frame_equal(output_df, expected_df)