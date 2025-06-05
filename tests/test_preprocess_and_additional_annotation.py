import pandas as pd
from src.preprocess_and_additional_annotation import (
    drop_abnormal_mapped_transcripts,
    cording_information_annotator,
    flame_information_annotator,
    max_min_exon_count_annotator,
    variant_count_annotator
)

def test_drop_abnormal_mapped_transcripts():
    # _GL や_MUの異常な染色体にマッピングされたトランスクリプトを削除
    input_data = pd.DataFrame({
    "geneName": ["gene1", "gene1", "gene2","gene2", "gene3","gene3", "gene4", "gene5"],
    "chrom": ["chr1", "chr1", "chr2", "chr2","chr3","chr3", "chrX_GL", "chr1_MU"],
    })
    output_data = drop_abnormal_mapped_transcripts(input_data)
    expected_output = pd.DataFrame({
    "geneName": ["gene1", "gene1", "gene2", "gene2", "gene3", "gene3"],
    "chrom": ["chr1", "chr1", "chr2", "chr2", "chr3", "chr3"],
    })
    pd.testing.assert_frame_equal(output_data, expected_output)

def test_cording_information_annotator():
    #遺伝子がcordingかnon-codingかを判定する
    input_data = pd.DataFrame({
    "cdsStart": [0, 100, 200, 500],
    "cdsEnd": [0, 300, 400, 500],
    })
    expected_output = pd.DataFrame({
        "cdsStart": [0, 100, 200, 500],
        "cdsEnd": [0, 300, 400, 500],
        "coding": ["non-coding", "coding", "coding", "non-coding"]
    })
    expected_output["coding"] = expected_output["coding"].astype("category")
    output_data = cording_information_annotator(input_data)
    pd.testing.assert_frame_equal(output_data, expected_output)

def test_flame_information_annotator():
    # exonの長さに基づいて、in-flame(mod3=0)かout-flame(mod3=1.2)かを判定する
    input_data = pd.DataFrame({
        "exonlengths": [[120, 200], [150, 250, 300]]
    })
    expected_output = pd.DataFrame({
        "exonlengths": [[120, 200], [150, 250, 300]],
        "flame": [
            ["in-flame", "out-flame"],
            ["in-flame", "out-flame", "in-flame"]
        ]
    })
    output_data = flame_information_annotator(input_data)
    pd.testing.assert_frame_equal(output_data, expected_output)

def test_max_min_exon_count_annotator():
    # 各遺伝子の最大および最小のエクソン数を追加する
    input_data = pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene2", "gene2", "gene3"],
        "exonCount": [2, 3, 1, 4, 2]
    })
    expected_output = pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene2", "gene2", "gene3"],
        "exonCount": [2, 3, 1, 4, 2],
        "max_exon_count": [3, 3, 4, 4, 2],
        "min_exon_count": [2, 2, 1, 1, 2]
    })
    output_data = max_min_exon_count_annotator(input_data)
    pd.testing.assert_frame_equal(output_data, expected_output)

def test_variant_count_annotator():
    # 各遺伝子のバリアント数を追加する
    input_data = pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene2", "gene2", "gene3"],
        "name": ["1variant1", "1variant2", "2variant1", "2variant2", "3variant1"]
    })
    expected_output = pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene2", "gene2", "gene3"],
        "name": ["1variant1", "1variant2", "2variant1", "2variant2", "3variant1"],
        "variant_count": [2,2,2,2,1],
    })
    output_data = variant_count_annotator(input_data)

    pd.testing.assert_frame_equal(output_data, expected_output)

