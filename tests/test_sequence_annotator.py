import pandas as pd

from altex_aid.sequence_annotator import (
    annotate_sequence_to_bed,
    join_sequence_to_single_exon_df,
)

TEST_FASTA_PATH = "tests/data/test.fa"
"""
テスト用のFASTAファイルの内容:
>chr1
AGCTAGCTAGCTAGCTAGCT
>chr2
GATTACAGATTACAGATTACA
つまり、chr1の0,4の範囲は"AGCT"、chr2の5,12の範囲は"CAGATTA"となる
"""


def test_annotate_sequence_to_bed():
    """
    FASTAファイルから正しく配列を抽出し、DataFrameに追加できるかテストする。
    プラス鎖とマイナス鎖の両方をテストする。
    """
    input_data = pd.read_csv("tests/data/test_bed_input.tsv", sep="\t", header=None)
    """
    テスト用のBED形式のデータ:
    chr1	0	4	gene1	0	+
    chr2	5	12	gene2	1	-
    """
    expected_output = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "chromStart": [0, 5],
            "chromEnd": [4, 12],
            "name": ["gene1", "gene2"],
            "score": [0, 1],
            "strand": ["+", "-"],
            "sequence": ["AGCT", "TAATCTG"],
            # chr1の0-4は"AGCT"、chr2の5-12は"CAGATTA"で、2列目はstrandが"-"なので相補鎖の"TAATCTG"が出力される
        }
    )
    output_data = annotate_sequence_to_bed(input_data, TEST_FASTA_PATH)
    print(output_data)
    pd.testing.assert_frame_equal(output_data, expected_output)


def test_join_sequence_to_single_exon_df():
    """
    acceptor_sequenceまたはdonor_sequenceをsingle_exon_dfに正しく追加できるかテストする。
    """

    # テスト用のsingle_exon_df
    single_exon_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "chromStart": [0, 5],
            "chromEnd": [4, 12],
            "name": ["gene1", "gene2"],
            "score": [0, 1],
            "strand": ["+", "-"],
        }
    )

    # --- acceptorの場合のテスト ---
    acceptor_bed_with_sequences = pd.DataFrame(
        {"score": [0, 1], "sequence": ["ATGC", "GTCTAAT"]}
    )
    donor_bed_with_sequences = pd.DataFrame(
        {"score": [0, 1], "sequence": ["TACG", "TAGATTA"]}
    )

    expected_output = single_exon_df.copy()
    expected_output = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "chromStart": [0, 5],
            "chromEnd": [4, 12],
            "name": ["gene1", "gene2"],
            "score": [0, 1],
            "strand": ["+", "-"],
            "acceptor_sequence": ["ATGC", "GTCTAAT"],
            "donor_sequence": ["TACG", "TAGATTA"],
        }
    )

    output_data = join_sequence_to_single_exon_df(
        single_exon_df,
        acceptor_bed_with_sequences,
        donor_bed_with_sequences,
    )
    pd.testing.assert_frame_equal(output_data, expected_output)
