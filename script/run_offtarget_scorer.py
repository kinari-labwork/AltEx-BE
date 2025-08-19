import pandas as pd
from pathlib import Path
from altex_aid.offtarget_scorer import (
    add_crisprdirect_url_to_df,
    calculate_offtarget_site_count_simple
)

exploded_sgrna_df = pd.read_pickle("data/formatted_exploded_sgrna_df.pkl")
assembly_name = "mm39"

fasta_path = Path("data/mm39/mm39.fa")

exploded_sgrna_df_with_crispr_direct = add_crisprdirect_url_to_df(exploded_sgrna_df, assembly_name)
exploded_sgrna_df_with_offtarget_counts = calculate_offtarget_site_count_simple(exploded_sgrna_df_with_crispr_direct, fasta_path)

exploded_sgrna_df_with_offtarget_counts.to_pickle("data/exploded_sgrna_df_with_offtarget_counts.pkl")