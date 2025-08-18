import pandas as pd
import pickle
from altex_aid.sgrna_designer import (
    design_sgrna_for_base_editors,
    make_preset_base_editors,
    design_sgrna_for_base_editors_dict
)

target_exons_df_with_acceptor_and_donor_sequence = pd.read_pickle("data/target_exons_with_acceptor_and_donor_sequence.pkl")

base_editors = make_preset_base_editors()
base_editors_list = list(base_editors.values())
target_exons_df_with_sgrna = design_sgrna_for_base_editors(
    target_exon_df=target_exons_df_with_acceptor_and_donor_sequence,
    base_editors=base_editors_list
)

target_exons_df_with_sgrna_dict = design_sgrna_for_base_editors_dict(
    target_exon_df=target_exons_df_with_acceptor_and_donor_sequence,
    base_editors=base_editors_list
)

print(target_exons_df_with_sgrna.columns)

target_exons_df_with_sgrna.to_pickle("data/target_exons_with_sgrna.pkl")

with open("data/target_exons_with_sgrna_dict.pkl", "wb") as f:
    pickle.dump(target_exons_df_with_sgrna_dict, f)
