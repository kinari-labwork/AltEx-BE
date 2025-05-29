import pandas as pd

from cassette_exon_extraction import classify_exons_per_gene


def test_classify_exons_per_gene():
        input_data = pd.DataFrame({
            "geneName": ["gene1", "gene1", "gene1","gene2","gene2"],
            "name": ["1-transcript1", "1-transcript2", "1-transcript3","2-transcript1", "2-transcript2"],
            "exons":[
                [(0, 100), (150, 200), (250, 300),],                  
                [(0, 100), (250, 300)],                               # (150,200) がスキップ
                [(0, 100), (110, 120), (150,200), (250, 300)],        # (110,120) のユニークなエキソンを持っている (ここまでが gene1)
                [(100,200), (250, 300)],                              # gene2のtranscript1
                [(100,200), (250, 350)],                             # gene2のtranscript2 (250,350のa3ssが生じている)
                ],   
            })
        expected_output = pd.DataFrame({
            "geneName": ["gene1", "gene1", "gene1", "gene2", "gene2"],
            "name": ["1-transcript1", "1-transcript2", "1-transcript3","2-transcript1", "2-transcript2"],
            "exons":[  
                [(0, 100), (150, 200), (250, 300),],                 
                [(0, 100), (250, 300)],                               # (150,200) がスキップ
                [(0, 100), (110, 120), (150,200), (250, 300)],        # (110,120) のユニークなエキソンを持っている (ここまでが gene1)
                [(100,200), (250, 300)],                              # gene2のtranscript1
                [(100,200), (250, 350)],                             # gene2のtranscript2 (250,350のa3ssが生じている)
                ],   
            "exontype": [
                ["constitutive", "cassette", "constitutive"],
                ["constitutive","constitutive"],
                ["constitutive", "unique", "cassette", "constitutive"],
                ["constitutive", "a3ss"],
                ["constitutive", "a3ss"],
                ],
            })
        output_data = classify_exons_per_gene(input_data)
        output_data = output_data[expected_output.columns]
        pd.testing.assert_frame_equal(output_data, expected_output)


