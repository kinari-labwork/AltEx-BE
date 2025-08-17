import pandas as pd
import numpy as np
from pathlib import Path
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
    テスト用のFastaファイル
    >chr1_test
    CCCCCGGGGATTACAGATTACAGATTACAGGGGGGATTACAGATTACAGATTACACCCCC
    >chr2_test
    AAAAAGCTAGCTAGCTAGCTAGCTATTTTTGGGGATTACAGATTACAGATTACA
    """

    fasta_path = Path("tests/data/test2.fa")

    input_df = pd.DataFrame({
        "uuid": ["id_A", "id_B", "id_C", "id_D"],
        "sgrna_target_sequence": [
            "GGGGATTACAGATTACAGATTACA",      # 23-mer, 3回出現
            "AAAAAAAAAAAAAAAAAAAAAAAA",      # 23-mer, 0回出現
            "GGGGATTACAGATTACAGATTACA",      # 重複ケース
            pd.NA
        ]
    })
    
    expected_df = pd.DataFrame({
        "uuid": ["id_A", "id_B", "id_C", "id_D"],
        "sgrna_target_sequence": [
            "GGGGATTACAGATTACAGATTACA",
            "AAAAAAAAAAAAAAAAAAAAAAAA",
            "GGGGATTACAGATTACAGATTACA",
            pd.NA
        ],
        "pam+20bp_exact_match_count": [3.0, 0.0, 3.0, np.nan]
    })

    output_df = calculate_offtarget_site_count(input_df, fasta_path)

    # Assert: 結果を検証
    pd.testing.assert_frame_equal(output_df, expected_df)