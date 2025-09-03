import pandas as pd

from altex_be.target_exon_extractor import (
    wrap_extract_target_exon
)

classified_refflat = pd.read_pickle("data/classified_exon_refflat.pkl")

splice_acceptor_single_exon_df, splice_donor_single_exon_df, exploded_classified_refflat = wrap_extract_target_exon(classified_refflat)

# BED形式で保存するためにタブ区切りで保存
splice_acceptor_single_exon_df.to_csv(
    "data/splice_acceptor_single_exon.bed", sep="\t", index=False, header=False
)
splice_donor_single_exon_df.to_csv(
    "data/splice_donor_single_exon.bed", sep="\t", index=False, header=False
)
exploded_classified_refflat.to_pickle("data/exploded_classified_refflat.pkl")
