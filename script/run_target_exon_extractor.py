import pandas as pd

from altex_aid.target_exon_extractor import (
    extract_splice_acceptor_regions,
    extract_splice_donor_regions,
    extract_target_exon,
)

classified_refflat = pd.read_pickle("data/classified_exon_refflat.pkl")

# BED形式も0base-start, 1base-endであるため、refFlatのexonStartsとexonEndsをそのまま使用する
target_exon_df = extract_target_exon(classified_refflat)
# SAから+-25bpの範囲を指定したdfの作成
splice_acceptor_single_exon_df = extract_splice_acceptor_regions(target_exon_df, 25)
# SDから+-25bpの範囲を指定したdfの作成
splice_donor_single_exon_df = extract_splice_donor_regions(target_exon_df, 25)

# BED形式で保存するためにタブ区切りで保存
target_exon_df.to_pickle("data/skipped_or_unique_exon_df.pkl")
splice_acceptor_single_exon_df.to_csv(
    "data/splice_acceptor_single_exon.bed", sep="\t", index=False, header=False
)
splice_donor_single_exon_df.to_csv(
    "data/splice_donor_single_exon.bed", sep="\t", index=False, header=False
)
