import pandas as pd

from cassette_exon_extraction import modify_refFlat

def test_modify_refFlat():
        input_data = pd.DataFrame({
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
        "exonEnds": ["90,200,300", "150,250,350", "50,150"]
    })
        expected_output = pd.DataFrame({
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
            [(0, 50), (50, 150)]
            ],
        "exonlengths": [
            [90, 100, 100],
            [150, 100, 100],
            [50, 100]
        ]   
    })
        output_data = modify_refFlat(input_data)
        output_data = output_data[expected_output.columns]
        pd.testing.assert_frame_equal(output_data, expected_output)

