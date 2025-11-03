import pandas as pd
import pickle
import time
from altex_be.sgrna_designer import (
    BaseEditor,
    design_sgrna_for_base_editors,
    design_sgrna_for_base_editors_dict,
)

target_exons_df_with_acceptor_and_donor_sequence = pd.read_pickle("data/target_exons_with_acceptor_and_donor_sequence.pkl")

base_editors = {}
base_editors["be4max-ng"] = (
    BaseEditor(
        base_editor_name="be4max-ng",
        pam_sequence="ng",
        editing_window_start_in_grna= 12,
        editing_window_end_in_grna= 17,
        base_editor_type= "cbe"
    )
)
base_editors["target-aid-ng"] = (
    BaseEditor(
        base_editor_name="target-aid-ng",
        pam_sequence="ng",
        editing_window_start_in_grna=17,
        editing_window_end_in_grna=19,
        base_editor_type= "cbe"
)
)
base_editors["abe8e-ng"] = (
    BaseEditor(
        base_editor_name="abe8e-ng",
        pam_sequence="ng",
        editing_window_start_in_grna=12,
        editing_window_end_in_grna=17,
        base_editor_type= "abe"
    )
)

t1 = time.time()

target_exons_df_with_sgrna = design_sgrna_for_base_editors(
    target_exon_df=target_exons_df_with_acceptor_and_donor_sequence,
    base_editors=base_editors
)

t2 = time.time()
print(f"Design sgRNA for base editors took {t2 - t1} seconds.")

t1 = time.time()
target_exons_df_with_sgrna_dict = design_sgrna_for_base_editors_dict(
    target_exon_df=target_exons_df_with_acceptor_and_donor_sequence,
    base_editors=base_editors
)
print(target_exons_df_with_sgrna.columns)
t2 = time.time()
print(f"Design sgRNA for base editors (dict version) took {t2 - t1} seconds.")

target_exons_df_with_sgrna.to_pickle("data/target_exons_with_sgrna.pkl")

with open("data/target_exons_with_sgrna_dict.pkl", "wb") as f:
    pickle.dump(target_exons_df_with_sgrna_dict, f)