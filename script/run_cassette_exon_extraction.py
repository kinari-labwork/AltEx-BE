import pandas as pd
import os
import sys
P=print

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))
from cassette_exon_extraction import modify_refFlat, classify_exons_per_gene, flip_a3ss_a5ss_in_minus_strand

refflat = pd.read_csv(
    "data/refFlat.txt",
    sep="\t",
    header=None,
    names=[
        "geneName",
        "name",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds"])

# "exonStarts" と　"exonEnds"を扱いやすいように (start, end) のリストに変換する
refflat_modify = modify_refFlat(refflat)
refflat_modify.to_csv("data/refFlat_modify.csv", index=False)
P(refflat_modify.head())

# 各遺伝子ごとにエクソンの分類を行い、"exontype"列に追加する
exon_classification = classify_exons_per_gene(refflat_modify)

# refflatは+ strandベースでexonのstart endが決まっているので、strand = - の遺伝子のa5ssとa3ssを入れ替える
exon_classification = flip_a3ss_a5ss_in_minus_strand(exon_classification)

#データフレームを保存する
exon_classification.to_csv("data/exon_classification.csv", index=False)
