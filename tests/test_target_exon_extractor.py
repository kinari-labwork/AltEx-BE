import pandas as pd

from altex_be.target_exon_extractor import (
    extract_splice_acceptor_regions,
    extract_splice_donor_regions,
    explode_classified_refflat,
    format_classified_refflat_to_bed,
)


def test_extract_target_exon():
    input_data = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2", "chr3"],
            "strand": ["+", "+", "-", "-", "+"],
            "exonStarts": [[100, 200, 300], [100, 400], [500, 600], [700, 850], [900]],
            "exonEnds": [[150, 250, 350], [150, 450], [550, 650], [750, 850], [1000]],
            "exons": [
                ["exon1", "exon2", "exon3"],
                ["exon1", "exon2"],
                ["exon1", "exon2"],
                ["exon1", "exon2"],
                ["exon1"],
            ],
            "exontype": [
                ["skipped", "constitutive", "skipped"],
                ["skipped", "constitutive"],
                ["a3ss-long", "constitutive"],
                ["a5ss-long", "constitutive"],
                ["constitutive"],
            ],
            "exonlengths": [
                [50, 100, 50],
                [50, 50],
                [100, 100],
                [100, 100],
                [1000],
            ],
            "flame":[
                ["in-flame", "in-flame", "out-flame"],
                ["in-flame", "out-flame"],
                ["in-flame", "out-flame"],
                ["in-flame", "out-flame"],
                ["in-flame"],
            ],
            "cds_info": [
                ["cds_start", "cds_end", "cds_start"],
                ["cds_start", "cds_end"],
                ["cds_start", "cds_end"],
                ["cds_start", "cds_end"],
                ["cds_start"],
            ],
            "exon_position": [
                ["first", "internal", "last"],
                ["first", "last"],
                ["first", "last"],
                ["first", "last"],
                ["single"],
            ],
        }
    )

    output_data = explode_classified_refflat(input_data, target_exon="alternative_exons")
    output_data = format_classified_refflat_to_bed(output_data)
    expected_output = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "exonStarts": [100, 300, 500, 700],
            "exonEnds": [150, 350, 550, 750],
            "score": [0, 0, 0, 0],
            # ここでは、ダミーのUUIDを使用（テストではoutput, expected_outputからは除外して比較する）
            "name": ["UUID1", "UUID2", "UUID3", "UUID4"],
            "strand": ["+", "+", "-", "-"],
            "exontype": ["skipped", "skipped", "a3ss-long", "a5ss-long"],
            "exon_position": ["first", "last", "first", "first"],
        }
    )
    print(output_data)
    print(expected_output)
    pd.testing.assert_frame_equal(
        output_data.drop(columns=["name"]).reset_index(drop=True), #UUIDは毎回異なる値が生成されるため、比較から除外
        expected_output.drop(columns=["name"]).reset_index(drop=True)
    )


def test_extract_splice_acceptor_regions():
    input_data = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "exonStarts": [100, 500],
            "exonEnds": [150, 550],
            "name" : ["UUID1", "UUID2"],
            "score": [0, 0],
            "strand": ["+", "-"],
            "exontype": ["skipped", "unique"],
            "exon_position": ["first", "last"],
        }
    )
    window = 25
    output_data = extract_splice_acceptor_regions(input_data, window)
    expected_output = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "chromStart": [75, 525],
            "chromEnd": [125, 575],
            "name": ["UUID1", "UUID2"],
            "score": [0, 0],
            "strand": ["+", "-"],
        }
    )
    pd.testing.assert_frame_equal(
        output_data.reset_index(drop=True), expected_output.reset_index(drop=True)
    )


def test_extract_splice_donor_regions():
    input_data = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "exonStarts": [100, 500],
            "exonEnds": [150, 550],
            "name" : ["UUID1", "UUID2"],
            "score": [0, 0],
            "strand": ["+", "-"],
            "exontype": ["skipped", "unique"],
            "exon_position": ["first", "last"],
        }
    )
    window = 25
    output_data = extract_splice_donor_regions(input_data, window)
    expected_output = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "chromStart": [125, 475],
            "chromEnd": [175, 525],
            "name": ["UUID1", "UUID2"],
            "score": [0, 0],
            "strand": ["+", "-"],
        }
    )
    pd.testing.assert_frame_equal(
        output_data.reset_index(drop=True), expected_output.reset_index(drop=True)
    )