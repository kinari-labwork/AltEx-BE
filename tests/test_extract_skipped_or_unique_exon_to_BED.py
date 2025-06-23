import pandas as pd
from extract_skipped_or_unique_exon_to_BED import extract_skipped_or_unique_exon, format_to_single_exon_bed, extract_splice_acceptor_regions, extract_splice_donor_regions

def test_extract_skipped_or_unique_exon():
    input_data = pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene2", "gene2","gene3"],
        "name": ["1transcript1", "1transcript2", "2transcript1", "2transcript2","3transcript1"],
        "chrom": ["chr1", "chr1", "chr2", "chr2","chr3"],
        "strand": ["+", "+", "-", "-","+"],
        "exonStarts": [[100, 200, 300], [100, 400], [500, 600], [700,850],[900]],
        "exonEnds": [[150, 250, 350], [150, 450], [550, 650], [750, 850],[1000]],
        "exontype": [["skipped","constitutive","skipped"], ["skipped","constitutive"], ["skipped","constitutive"], ["unique","constitutive"],["constitutive"]],
        "exon_position":[["first","internal","last"],["first","last"],["first","last"],["first","last"],["single"]]
    })

    output_data = extract_skipped_or_unique_exon(input_data)
    expected_output = pd.DataFrame({
        "geneName": ["gene1","gene1","gene2", "gene2"],
        "name": ["1transcript1", "1transcript1","2transcript1", "2transcript2"],
        "chrom": ["chr1", "chr1", "chr2", "chr2"],
        "strand": ["+","+", "-", "-"],
        "exonStarts": [100, 300, 500, 700],
        "exonEnds": [150, 350, 550, 750],
        "exontype": ["skipped","skipped", "skipped","unique"],
        "exon_position": ["first","last","first", "first"],
        "score": [0, 1, 2, 3]
    })
    print(output_data)
    print(expected_output)
    pd.testing.assert_frame_equal(output_data, expected_output)

def test_extract_splice_acceptor_regions():
    input_data = pd.DataFrame({
        "geneName": ["gene1", "gene2"],
        "chrom": ["chr1", "chr2"],
        "strand": ["+", "-"],
        "exonStarts": [100, 500],
        "exonEnds": [150, 550],
        "score": [0, 1]
    })
    window = 25
    output_data = extract_splice_acceptor_regions(input_data, window)
    expected_output = pd.DataFrame({
        "geneName": ["gene1", "gene2"],
        "chrom": ["chr1", "chr2"],
        "strand": ["+", "-"],
        "chromStart": [75, 525],
        "chromEnd": [125, 575],
        "score": [0, 1]
    })
    pd.testing.assert_frame_equal(output_data.reset_index(drop=True), expected_output.reset_index(drop=True))

def test_extract_splice_donor_regions():
    input_data = pd.DataFrame({
        "geneName": ["gene1", "gene2"],
        "chrom": ["chr1", "chr2"],
        "strand": ["+", "-"],
        "exonStarts": [100, 500],
        "exonEnds": [150, 550],
        "score": [0, 1]
    })
    window = 25
    output_data = extract_splice_donor_regions(input_data, window)
    expected_output = pd.DataFrame({
        "geneName": ["gene1", "gene2"],
        "chrom": ["chr1", "chr2"],
        "strand": ["+", "-"],
        "chromStart": [125, 475],
        "chromEnd": [175, 525],
        "score": [0, 1]
    })
    pd.testing.assert_frame_equal(output_data.reset_index(drop=True), expected_output.reset_index(drop=True))

def test_format_to_single_exon_bed():
    input_data = pd.DataFrame({
        "geneName": ["gene1","gene2", "gene2"],
        "chrom": ["chr1", "chr2", "chr2"],
        "strand": ["+", "-", "-"],
        "chromStart": [100, 500, 700],
        "chromEnd": [150, 550, 750],
        "exontype": ["skipped","skipped","unique"],
        "exon_position": ["first", "last", "first"],
        "score" : [0, 1, 2]
    })
    output_data = format_to_single_exon_bed(input_data)
    expected_output = pd.DataFrame({
        "chrom": ["chr1", "chr2", "chr2"],
        "chromStart": [100, 500, 700],
        "chromEnd": [150, 550, 750],
        "name": ["gene1", "gene2", "gene2"],
        "score": [0, 1, 2],
        "strand": ["+", "-", "-"]
    })
    pd.testing.assert_frame_equal(output_data.reset_index(drop=True), expected_output.reset_index(drop=True))