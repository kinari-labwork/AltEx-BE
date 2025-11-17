import pandas as pd

from altex_be.refflat_preprocessor import (
    add_exon_position_flags,
    annotate_cording_information,
    annotate_flame_information,
    annotate_variant_count,
    calculate_exon_lengths,
    drop_abnormal_mapped_transcripts,
    parse_exon_coordinates,
)


def test_parse_exon_coordinates():
    input_data = pd.DataFrame(
        {
            "geneName": ["gene1", "gene1", "gene2"],
            "name": ["1-transcript1", "1-transcript2", "2-transcript1"],
            "chrom": ["chr1", "chr1", "chr2"],
            "strand": ["+", "+", "-"],
            "txStart": [0, 0, 0],
            "txEnd": [1000, 1000, 500],
            "cdsStart": [0, 0, 0],
            "cdsEnd": [1000, 1000, 500],
            "exonCount": [3, 3, 2],
            "exonStarts": ["0,100,200", "0,150,250", "0,50"],
            "exonEnds": ["90,200,300", "150,250,350", "50,150"],
        }
    )
    expected_output = pd.DataFrame(
        {
            "geneName": ["gene1", "gene1", "gene2"],
            "name": ["1-transcript1", "1-transcript2", "2-transcript1"],
            "chrom": ["chr1", "chr1", "chr2"],
            "strand": ["+", "+", "-"],
            "txStart": [0, 0, 0],
            "txEnd": [1000, 1000, 500],
            "cdsStart": [0, 0, 0],
            "cdsEnd": [1000, 1000, 500],
            "exonCount": [3, 3, 2],
            "exons": [
                [(0, 90), (100, 200), (200, 300)],
                [(0, 150), (150, 250), (250, 350)],
                [(0, 50), (50, 150)],
            ],
        }
    )
    output_data = parse_exon_coordinates(input_data)
    output_data = output_data[expected_output.columns]
    pd.testing.assert_frame_equal(output_data, expected_output)


def test_calculate_exon_lengths():
    # exonの長さを計算する
    input_data = pd.DataFrame(
        {
            "geneName": ["gene1", "gene2"],
            "exonStarts": [[0, 100, 150], [50, 200]],
            "exonEnds": [[100, 200, 300], [100, 300]],
        }
    )
    expected_output = pd.DataFrame(
        {
            "geneName": ["gene1", "gene2"],
            "exonStarts": [[0, 100, 150], [50, 200]],
            "exonEnds": [[100, 200, 300], [100, 300]],
            "exonlengths": [[100, 100, 150], [50, 100]],
        }
    )
    output_data = calculate_exon_lengths(input_data)
    pd.testing.assert_frame_equal(output_data, expected_output)


def test_drop_abnormal_mapped_transcripts():
    # _GL や_MUの異常な染色体にマッピングされたトランスクリプトを削除
    input_data = pd.DataFrame(
        {
            "geneName": [
                "gene1",
                "gene1",
                "gene2",
                "gene2",
                "gene3",
                "gene3",
                "gene4",
                "gene5",
            ],
            "chrom": [
                "chr1",
                "chr1",
                "chr2",
                "chr2",
                "chr3",
                "chr3",
                "chrX_random",
                "chr1_alt",
            ],
        }
    )
    output_data = drop_abnormal_mapped_transcripts(input_data)
    expected_output = pd.DataFrame(
        {
            "geneName": ["gene1", "gene1", "gene2", "gene2", "gene3", "gene3"],
            "chrom": ["chr1", "chr1", "chr2", "chr2", "chr3", "chr3"],
        }
    )
    pd.testing.assert_frame_equal(output_data, expected_output)


def test_cording_information_annotator():
    # 遺伝子がcordingかnon-codingかを判定する
    input_data = pd.DataFrame(
        {
            "name": ["NR_001", "NM_001", "NM_002", "NR_002"],
            "cdsStart": [0, 100, 50, 0],
            "cdsEnd": [0, 200, 450, 0],
        }
    )
    expected_output = pd.DataFrame(
        {
            "name": ["NR_001", "NM_001", "NM_002", "NR_002"],
            "cdsStart": [0, 100, 50, 0],
            "cdsEnd": [0, 200, 450, 0],
            "coding": ["non-coding", "coding", "coding", "non-coding"],
        }
    )
    expected_output["coding"] = expected_output["coding"].astype("category")
    output_data_refflat = annotate_cording_information(input_data, gtf_flag=False)
    pd.testing.assert_frame_equal(output_data_refflat, expected_output)
    output_data_gtf = annotate_cording_information(input_data, gtf_flag=True)
    pd.testing.assert_frame_equal(output_data_gtf, expected_output)


def test_flame_information_annotator():
    # exonの長さに基づいて、in-flame(mod3=0)かout-flame(mod3=1.2)かを判定する
    input_data = pd.DataFrame({"exonlengths": [[120, 200], [150, 250, 300]]})
    expected_output = pd.DataFrame(
        {
            "exonlengths": [[120, 200], [150, 250, 300]],
            "flame": [["in-flame", "out-flame"], ["in-flame", "out-flame", "in-flame"]],
        }
    )
    output_data = annotate_flame_information(input_data)
    pd.testing.assert_frame_equal(output_data, expected_output)


def test_variant_count_annotator():
    # 各遺伝子のバリアント数を追加する
    input_data = pd.DataFrame(
        {
            "geneName": ["gene1", "gene1", "gene2", "gene2", "gene3"],
            "name": ["1variant1", "1variant2", "2variant1", "2variant2", "3variant1"],
        }
    )
    expected_output = pd.DataFrame(
        {
            "geneName": ["gene1", "gene1", "gene2", "gene2", "gene3"],
            "name": ["1variant1", "1variant2", "2variant1", "2variant2", "3variant1"],
            "variant_count": [2, 2, 2, 2, 1],
        }
    )
    output_data = annotate_variant_count(input_data)
    pd.testing.assert_frame_equal(output_data, expected_output)


def test_add_exon_position_flags():
    input_data = pd.DataFrame(
        {
            "geneName": ["gene1", "gene2", "gene3"],
            "exonStarts": [[100, 200, 300], [400, 500], [200]],
            "strand": ["+", "-", "+"],
        }
    )
    expected_output = pd.DataFrame(
        {
            "geneName": ["gene1", "gene2", "gene3"],
            "exonStarts": [[100, 200, 300], [400, 500], [200]],
            "strand": ["+", "-", "+"],
            "exon_position": [
                ["first", "internal", "last"],
                ["last", "first"],
                ["single"],
            ],
        }
    )
    output_data = add_exon_position_flags(input_data)
    pd.testing.assert_frame_equal(output_data, expected_output)
