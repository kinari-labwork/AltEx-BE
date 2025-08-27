import pandas as pd

from altex_be.splicing_event_classifier import (
    classify_splicing_events_per_gene,
    flip_a3ss_a5ss_on_minus_strand,
)

refflat = pd.read_pickle("data/processed_refflat.pkl")

# 各遺伝子ごとにエクソンの分類を行い、"exontype"列に追加する
classified_refflat = classify_splicing_events_per_gene(refflat)

# refflatは+ strandベースでexonのstart endが決まっているので、strand = - の遺伝子のa5ssとa3ssを入れ替える
classified_refflat = flip_a3ss_a5ss_on_minus_strand(classified_refflat)

# データフレームを保存する
classified_refflat.to_pickle("data/classified_exon_refflat.pkl")
