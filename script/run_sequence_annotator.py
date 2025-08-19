import pandas as pd

from altex_aid.sequence_annotator import (
    annotate_sequence_to_bed,
    join_sequence_to_single_exon_df,
)

target_exons_df = pd.read_pickle("data/exploded_classified_refflat.pkl")
splice_acceptor_bed = pd.read_csv(
    "data/splice_acceptor_single_exon.bed", sep="\t", header=None
)
splice_donor_bed = pd.read_csv(
    "data/splice_donor_single_exon.bed", sep="\t", header=None
)
annotation_genome = "mm39"
fasta_path = f"data/{annotation_genome}/{annotation_genome}.fa"

# 塩基配列を取得して、BED形式のデータに追加
splice_acceptor_single_exon_df = annotate_sequence_to_bed(
    bed=splice_acceptor_bed, fasta_path=fasta_path
)
splice_donor_single_exon_df = annotate_sequence_to_bed(
    bed=splice_donor_bed, fasta_path=fasta_path
)

# 取得した塩基配列をtarget_exon_dfに結合
target_exon_df_with_acceptor_and_donor_sequence = join_sequence_to_single_exon_df(
    single_exon_df=target_exons_df,
    acceptor_bed_with_sequences=splice_acceptor_single_exon_df,
    donor_bed_with_sequences=splice_donor_single_exon_df,
)
# 結果を保存
target_exon_df_with_acceptor_and_donor_sequence.to_pickle(
    "data/target_exons_with_acceptor_and_donor_sequence.pkl"
)
