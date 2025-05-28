import sys
import unittest
import pandas as pd
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

from cassette_exon_extraction import classify_exon_type

class TestCassetteExonClassification(unittest.TestCase):
    def setUp(self):
        # 準備：複数のtranscriptを用意（3トランスクリプト）
        self.transcripts = [
            [(0, 100), (150, 200), (250, 300)],       # 全部あり（完全型）
            [(0, 100), (250, 300)],                   # (150,200) がスキップ
            [(0, 100), (150, 200), (250, 300)],        # (150,200) ではなく5'のalteration
            [(0, 100), (180, 220), (250, 300)]        # (180,220) のユニークなエキソンを持っている
        ]

    def test_constitutive(self):
        # (0, 100) はすべてのtranscriptに含まれる
        label = classify_exon_type((0, 100), self.transcripts)
        self.assertEqual(label, "constitutive")

    def test_cassette_exact(self):
        # (150, 200) は transcript1にだけ存在、他では flanking されている
        label = classify_exon_type((150, 200), self.transcripts)
        self.assertEqual(label, "cassette")

    def test_alt_5ss(self):
        # (155, 200) が存在し、startがずれているため alt_5ss を拾いたい
        # NOTE: (150,200) に対してズレたものがいるので、片側一致
        label = classify_exon_type((155, 200), self.transcripts)
        self.assertEqual(label, "a3ss")

    def test_unique(self):
        # (180, 220) はどこにも含まれないが、ここでuniqueの条件にあわせるなら
        # 自分だけ含まれている場合（例外的にmanualで渡しても可）
        label = classify_exon_type((180, 220), self.transcripts)
        self.assertEqual(label, "unique")

    def test_fuzzy(self):
        # (151,199) は他のexonに±3以内で近い（たとえば fuzzy tolerance内）
        label = classify_exon_type((151, 199), self.transcripts, fuzzy_tolerance=3)
        self.assertEqual(label, "fuzzy")

if __name__ == "__main__":
    unittest.main()