import pandas as pd

from altex_aid.splice_event_classifier import classify_exon_type

class TestClassifyExonType:
    def setup_method(self):
        # 準備：複数のtranscriptを用意（3トランスクリプト）
        self.transcripts = [
            [(0, 100), (150, 200), (250, 300)],       # 全部あり（完全型）
            [(0, 100), (250, 300)],                   # (150,200) がスキップ
            [(0, 100), (150, 200), (250, 300)],        # (150,200) ではなく5'のalteration
            [(0, 100), (110,120), (250, 300)],
            [(0, 100), (250,300), (500, 600), (600, 700)]
        ]

    def test_constitutive(self):
        # (0, 100) はすべてのtranscriptに含まれる
        label = classify_exon_type((0, 100), self.transcripts)
        assert label == "constitutive"

    def test_cassette_exact(self):
        # (150, 200) は transcript1にだけ存在、他では flanking されている
        label = classify_exon_type((150, 200), self.transcripts)
        assert label == "skipped"

    def test_alt_3ss(self):
        # (155, 200) が存在し、startがずれているため a3ss を拾いたい
        label = classify_exon_type((155, 200), self.transcripts)
        assert label == "a3ss"     
    
    def test_alt_5ss(self):
        # (150,200) に対して(150,205)はendがずれているため a5ss を拾いたい
        label = classify_exon_type((150, 205), self.transcripts)
        assert label == "a5ss"
    def test_split(self):
        #(500,700)は(500,600)と(600,700)に分割されている
        label = classify_exon_type((500, 700), self.transcripts)
        assert label == "intron_retention"

    def test_unique(self):
        # (180, 220) はどこにも含まれない
        label = classify_exon_type((110, 120), self.transcripts)
        assert label == "unique"

    def test_overlap(self):
        # (151, 199)は(150,200)とオーバーラップしている 
        label = classify_exon_type((151, 199), self.transcripts)
        assert  label == "overlap"


from altex_aid.splice_event_classifier import classify_exons_per_gene

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
                [(100,200), (250, 350)],                             # gene2のtranscript2 (250,350のa5ssが生じている)
                ],   
            "exontype": [
                ["constitutive", "skipped", "constitutive"],
                ["constitutive","constitutive"],
                ["constitutive", "unique", "skipped", "constitutive"],
                ["constitutive", "a5ss"],
                ["constitutive", "a5ss"],
                ],
            })
        output_data = classify_exons_per_gene(input_data)
        output_data = output_data[expected_output.columns]
        pd.testing.assert_frame_equal(output_data, expected_output)


from altex_aid.splice_event_classifier import flip_a3ss_a5ss_in_minus_strand

def test_flip_a3ss_a5ss_in_minus_strand():
         input_data = pd.DataFrame({
            "strand": ["+", "+", "-", "-"],
            "exontype": [
                ["constitutive", "skipped", "a5ss"],
                ["constitutive", "overlap", "a3ss"],
                ["constitutive", "unique", "a5ss"],
                ["constitutive", "cassette", "a3ss"],
                ],
            }) 
         expected_output = pd.DataFrame({
              "strand": ["+", "+", "-", "-"],
              "exontype": [
                ["constitutive", "skipped", "a5ss"],  # plus strandは変化なし
                ["constitutive", "overlap", "a3ss"],
                ["constitutive", "unique", "a3ss"],  # minus strandでだけ、a3ssとa5ssが入れ替わっている
                ["constitutive", "cassette", "a5ss"],
                ],
            })
         output_data = flip_a3ss_a5ss_in_minus_strand(input_data)
         pd.testing.assert_frame_equal(output_data,expected_output)