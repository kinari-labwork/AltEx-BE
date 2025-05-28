import sys
import unittest
import pandas as pd
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

from cassette_exon_extraction import classify_exons_per_gene

class TestClassifyExonsPerGene(unittest.TestCase):
    def test_classify_exons_per_gene(self):
        input_data = pd.DataFrame({
            "geneName": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "name": ["1-transcript1", "1-transcript2", "1-transcript3", "1-transcript4","2-transcript1", "2-transcript2"],
            "chrom": ["chr1", "chr1", "chr1", "chr1", "chr2", "chr2"],
            "strand": ["+", "+", "+", "+", "-", "-"],
            "txStart": [0, 0, 0, 0, 1000, 1000],
            "txEnd": [1000, 1000, 500, 500, 2000, 2000],
            "cdsStart": [0, 0, 0, 0, 1000, 1000],
            "cdsEnd": [1000, 1000, 500, 500, 2000, 2000],
            "exonCount": [3, 2, 3, 3, 2, 2],
            "exons":[       
                [(0, 100), (150, 200), (250, 300)],       # 全部あり（完全型）
                [(0, 100), (250, 300)],                   # (150,200) がスキップ
                [(0, 100), (150, 200), (250, 300)],       # (150,200) ではなく5'のalteration
                [(0, 100), (180, 220), (250, 300)],       # (180,220) のユニークなエキソンを持っている (ここまでが gene1)
                [(100,200), (250, 300)],                  # gene2のtranscript1
                [(100, 200), (250, 300)],                 # gene2のtranscript2
                ],   
            "exonlengths": [
                [100, 50, 50],
                [100, 50],
                [100, 50, 50],
                [100, 40, 50],
                [100, 50],
                [100, 50]
                ]
            })
        expected_output = pd.DataFrame({
            "geneName": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2"],
            "name": ["1-transcript1", "1-transcript2", "1-transcript3", "1-transcript4","2-transcript1", "2-transcript2"],
            "chrom": ["chr1", "chr1", "chr1", "chr1", "chr2", "chr2"],
            "strand": ["+", "+", "+", "+", "-", "-"],
            "txStart": [0, 0, 0, 0, 1000, 1000],
            "txEnd": [1000, 1000, 500, 500, 2000, 2000],
            "cdsStart": [0, 0, 0, 0, 1000, 1000],
            "cdsEnd": [1000, 1000, 500, 500, 2000, 2000],
            "exonCount": [3, 2, 3, 3, 2, 2],
            "exons":[        
                [(0, 100), (150, 200), (250, 300)],       # 全部あり（完全型）
                [(0, 100), (250, 300)],                   # (150,200) がスキップ
                [(0, 100), (150, 200), (250, 300)],       # (150,200) ではなく5'のalteration
                [(0, 100), (180, 220), (250, 300)],       # (180,220) のユニークなエキソンを持っている (ここまでが gene1)
                [(100,200), (250, 300)],                  # gene2のtranscript1
                [(100, 200), (250, 300)],                 # gene2のtranscript2
                ],   
            #つまり、cassetteは(150,200)であり、uniqueは(180,220)である。
            "exonlengths": [
                [100, 50, 50],
                [100, 50],
                [100, 50, 50],
                [100, 40, 50],
                [100, 50],
                [100, 50]
                ],
            "exontype": [
                ["constitutive", "cassette", "constitutive"],
                ["constitutive","constitutive"],
                ["constitutive", "cassette", "constitutive"],
                ["constitutive", "unique", "constitutive"],
                ["constitutive", "constitutive"],
                ["constitutive", "constitutive"],
                ],
            })
        output_data = classify_exons_per_gene(input_data)
        output_data = output_data[expected_output.columns]
        pd.testing.assert_frame_equal(output_data, expected_output)

if __name__ == "__main__":
    unittest.main()