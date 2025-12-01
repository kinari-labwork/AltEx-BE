import pandas as pd
from pathlib import Path
from altex_be.offtarget_scorer import (
    add_crisprdirect_url_to_df,
    calculate_offtarget_site_count_ahocorasick,
    add_reversed_complement_sgrna_column
)

def test_add_crisprdirect_url_to_df():
    # 1. Arrange: テストデータの準備
    input_df = pd.DataFrame({
        "sgrna_target_sequence": ["GAT+TACA", "ATA+TATA"],
        "base_editor_pam_sequence": ["NGG", "NGG"],
    })
    
    expected_df = pd.DataFrame({
        "sgrna_target_sequence": ["GAT+TACA", "ATA+TATA"],
        "base_editor_pam_sequence": ["NGG", "NGG"],
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
    AAAAAGCTAGCTAGCTAGCTAGCTATTTTTGGGGATTACAGATTATCGATACA
    """

    fasta_path = Path("tests/data/test2.fa")

    input_df = pd.DataFrame({
        "strand": ["+", "+", "+", "-"],
        "uuid": ["id_A", "id_B", "id_C", "id_D"],
        "sgrna_target_sequence": [
            "GGG+GATTACAGATTACAGATTAC",      # 20-mer, 2回出現, 12-merは3回出現
            "AAA+AAAAAAAAAAAAAAAAAAAA",      # 20-mer, 0回出現
            "ggg+gattacagattacagattac",      # 小文字にした場合のテストケース
            "GTAATCTGTAATCTGTAATC+CCC"       # 20-mer, 逆相補化したもの。 2回出現, 12-merは3回出現
        ]
    })
    
    expected_df = pd.DataFrame({
        "strand": ["+", "+", "+", "-"],
        "uuid": ["id_A", "id_B", "id_C", "id_D"],
        "sgrna_target_sequence": [
            "GGG+GATTACAGATTACAGATTAC",
            "AAA+AAAAAAAAAAAAAAAAAAAA",
            "ggg+gattacagattacagattac",
            "GTAATCTGTAATCTGTAATC+CCC" 
        ],
        "pam+20bp_exact_match_count": [2, 0, 2, 2],
        "pam+12bp_exact_match_count": [3, 0, 3, 3]

    })

    input_df = add_reversed_complement_sgrna_column(input_df)
    output_df = calculate_offtarget_site_count_ahocorasick(input_df, fasta_path)
    print(output_df)

    # Assert: 結果を検証
    pd.testing.assert_frame_equal(output_df, expected_df)