import pandas as pd

from altex_be.splicing_event_classifier import (
    classify_splicing_event,
    classify_splicing_events_per_gene,
    flip_a3ss_a5ss_on_minus_strand,
)


class TestClassifyExonType:
    def setup_method(self):
        # 準備：複数のtranscriptを用意（3トランスクリプト）
        self.transcripts = [
            [(0, 100), (150, 200), (250, 300)],  # 全部あり（完全型）
            [(0, 100), (250, 300)],  # (150,200) がスキップ
            [(0, 100), (150, 200), (250, 300)],  # (150,200) ではなく5'のalteration
            [(0, 100), (110, 120), (250, 300)],
            [(0, 100), (250, 300), (500, 600), (600, 700)],
        ]

    def test_constitutive(self):
        # (0, 100) はすべてのtranscriptに含まれる
        label = classify_splicing_event((0, 100), self.transcripts)
        assert label == "constitutive"

    def test_cassette_exact(self):
        # (150, 200) は transcript1にだけ存在、他では flanking されている
        label = classify_splicing_event((150, 200), self.transcripts)
        assert label == "skipped"

    def test_alt_3ss_short(self):
        # (150, 200) が存在し、(155,200)はstartがずれている。 しかし、end を共有する(150,200)のほうが長いため a3ss-short を拾いたい
        label = classify_splicing_event((155, 200), self.transcripts)
        assert label == "a3ss-short"
    
    def test_alt_3ss_long(self):
        # (150, 200) が存在し、(140,200)はstartがずれている。 これは end が一致する他のエキソンのstartよりも小さいため a3ss-long を拾いたい
        label = classify_splicing_event((140, 200), self.transcripts)
        assert label == "a3ss-long"

    def test_alt_5ss_short(self):
        # (150,200) に対して(150,190)はendがずれている。しかし、startを 共有する(150,200)のほうが長いため a5ss-short を拾いたい
        label = classify_splicing_event((150, 190), self.transcripts)
        assert label == "a5ss-short"

    def test_alt_5ss_long(self):
        # (150,200) に対して(150,210)はendがずれている。これは start が一致する他のエキソンのendよりも大きいため、a5ss-long を拾いたい
        label = classify_splicing_event((150, 210), self.transcripts)
        assert label == "a5ss-long"
    def test_split(self):
        # (500,700)は(500,600)と(600,700)に分割されている
        label = classify_splicing_event((500, 700), self.transcripts)
        assert label == "intron_retention"

    def test_unique(self):
        # (180, 220) はどこにも含まれない
        label = classify_splicing_event((110, 120), self.transcripts)
        assert label == "unique"

    def test_overlap(self):
        # (151, 199)は(150,200)とオーバーラップしている
        label = classify_splicing_event((151, 199), self.transcripts)
        assert label == "overlap"


def test_classify_exons_per_gene():
    input_data = pd.DataFrame(
        {
            "geneName": ["gene1", "gene1", "gene1", "gene2", "gene2"],
            "name": [
                "1-transcript1",
                "1-transcript2",
                "1-transcript3",
                "2-transcript1",
                "2-transcript2",
            ],
            "exons": [
                [
                    (0, 100),
                    (150, 200),
                    (250, 300),
                ],
                [(0, 100), (250, 300)],  # (150,200) がスキップ
                [
                    (0, 100),
                    (110, 120),
                    (150, 200),
                    (250, 300),
                ],  # (110,120) のユニークなエキソンを持っている (ここまでが gene1)
                [(100, 200), (250, 300)],  # gene2のtranscript1
                [
                    (100, 200),
                    (250, 350),
                ],  # gene2のtranscript2 (250,350のa3ssが生じている)
            ],
        }
    )
    expected_output = pd.DataFrame(
        {
            "geneName": ["gene1", "gene1", "gene1", "gene2", "gene2"],
            "name": [
                "1-transcript1",
                "1-transcript2",
                "1-transcript3",
                "2-transcript1",
                "2-transcript2",
            ],
            "exons": [
                [
                    (0, 100),
                    (150, 200),
                    (250, 300),
                ],
                [
                    (0, 100), 
                    (250, 300)
                ],  # (150,200) がスキップ
                [
                    (0, 100),
                    (110, 120),
                    (150, 200),
                    (250, 300),
                ],  # (110,120) のユニークなエキソンを持っている (ここまでが gene1)
                [
                    (100, 200),
                    (250, 300)
                    ],  # gene2のtranscript1
                [
                    (100, 200),
                    (250, 350),
                ],  # gene2のtranscript2 (250,350のa5ssが生じている)
            ],
            "exontype": [
                ["constitutive", "skipped", "constitutive"],
                ["constitutive", "constitutive"],
                ["constitutive", "unique", "skipped", "constitutive"],
                ["constitutive", "a5ss-short"],
                ["constitutive", "a5ss-long"],
            ],
        }
    )
    output_data = classify_splicing_events_per_gene(input_data)
    output_data = output_data[expected_output.columns]
    pd.testing.assert_frame_equal(output_data, expected_output)


def test_flip_a3ss_a5ss_in_minus_strand():
    input_data = pd.DataFrame(
        {
            "strand": ["+", "+", "-", "-"],
            "exontype": [
                ["constitutive", "skipped", "a5ss-long"],
                ["constitutive", "overlap", "a3ss-short"],
                ["constitutive", "unique", "a5ss-short"],
                ["constitutive", "cassette", "a3ss-long"],
            ],
        }
    )
    expected_output = pd.DataFrame(
        {
            "strand": ["+", "+", "-", "-"],
            "exontype": [
                ["constitutive", "skipped", "a5ss-long"],  # plus strandは変化なし
                ["constitutive", "overlap", "a3ss-short"],
                ["constitutive","unique","a3ss-short",],  # minus strandでだけ、a3ssとa5ssが入れ替わっている
                ["constitutive", "cassette", "a5ss-long"],
            ],
        }
    )
    output_data = flip_a3ss_a5ss_on_minus_strand(input_data)
    pd.testing.assert_frame_equal(output_data, expected_output)
