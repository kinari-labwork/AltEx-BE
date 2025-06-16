import pandas as pd
from refflat_to_single_exon_BED import refflat_to_single_exon_df, format_to_single_exon_bed

def test_refflat_to_single_exon():
    # テスト用のデータフレームを作成
    input_data = pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene2", "gene2"],
        "name": ["1transcript1", "1transcript2", "2transcript1", "2transcript2"],
        "chrom": ["chr1", "chr1", "chr2", "chr2"],
        "strand": ["+", "+", "-", "-"],
        "exonStarts": [[100, 200], [100, 400], [500, 600], [700, 800]],
        "exonEnds": [[150, 250], [150, 450], [550, 650], [750, 850]],
        "exontype": [["skipped","constitutive"], ["skipped","constitutive"], ["skipped","constitutive"], ["unique","constitutive"]]
    })

    output_data = refflat_to_single_exon_df(input_data)
    expected_output = pd.DataFrame({
        "geneName": ["gene1","gene2", "gene2"],
        "name": ["1transcript1", "2transcript1", "2transcript2"],
        "chrom": ["chr1", "chr2", "chr2"],
        "strand": ["+", "-", "-"],
        "exonStarts": [100, 500, 700],
        "exonEnds": [150, 550, 750],
        "exontype": ["skipped","skipped","unique"]
    })

    pd.testing.assert_frame_equal(output_data, expected_output)

def test_format_to_single_exon_bed():
    input_data = pd.DataFrame({
        "geneName": ["gene1","gene2", "gene2"],
        "chrom": ["chr1", "chr2", "chr2"],
        "strand": ["+", "-", "-"],
        "exonStarts": [100, 500, 700],
        "exonEnds": [150, 550, 750],
        "exontype": ["skipped","skipped","unique"]
    })
    output_data = format_to_single_exon_bed(input_data)
    expected_output = pd.DataFrame({
        "chrom": ["chr1", "chr2", "chr2"],
        "chromStart": [100, 500, 700],
        "chromEnd": [150, 550, 750],
        "name": ["gene1", "gene2", "gene2"],
        "strand": ["+", "-", "-"]
    })
    pd.testing.assert_frame_equal(output_data.reset_index(drop=True), expected_output.reset_index(drop=True))