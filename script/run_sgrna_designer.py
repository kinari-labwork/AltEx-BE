import pandas as pd
from altex_aid.sgrna_designer import (
    design_sgrna_for_base_editors,
    make_default_base_editors,
)

target_exons_df_with_acceptor_and_donor_sequence = pd.read_pickle("data/target_exons_with_acceptor_and_donor_sequence.pkl")

base_editors = make_default_base_editors()
target_exons_df_with_sgrna = design_sgrna_for_base_editors(
    target_exon_df=target_exons_df_with_acceptor_and_donor_sequence,
    base_editors=base_editors
)



target_exons_df_with_sgrna.to_pickle("data/target_exons_with_sgrna.pkl")

