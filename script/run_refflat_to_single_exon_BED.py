import pandas as pd
from extract_skipped_or_unique_exon_to_BED import extract_skipped_or_unique_exon, format_to_single_exon_bed

data = pd.read_pickle("data/exon_classification_with_additional_info.pkl")

# BED形式も0base-start, 1base-endであるため、refFlatのexonStartsとexonEndsをそのまま使用する
single_exon_df = extract_skipped_or_unique_exon(data)
single_exon_bed = format_to_single_exon_bed(single_exon_df)
# BED形式で保存するためにタブ区切りで保存
single_exon_df.to_pickle("data/skipped_or_unique_exon_df.pkl")
single_exon_bed.to_csv("data/skipped_or_unique_exon.bed", sep="\t", index=False, header=True)