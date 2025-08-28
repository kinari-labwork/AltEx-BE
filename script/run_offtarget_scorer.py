import pandas as pd
from pathlib import Path
from altex_be.offtarget_scorer import (
    score_offtargets
)

exploded_sgrna_df = pd.read_pickle("data/formatted_exploded_sgrna_df.pkl")
assembly_name = "mm39"

fasta_path = Path("data/mm39/mm39.fa")

exploded_sgrna_df_with_offtarget_counts = score_offtargets(exploded_sgrna_df, assembly_name, fasta_path)

exploded_sgrna_df_with_offtarget_counts.to_pickle("data/exploded_sgrna_df_with_offtarget_counts.pkl")