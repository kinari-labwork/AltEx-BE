import pandas as pd
from pathlib import Path
import pickle
from altex_be.output_formatter import format_output
from altex_be.sgrna_designer import BaseEditor, make_preset_base_editors

with open("data/target_exons_with_sgrna_dict.pkl", "rb") as f:
    target_exons_df_with_sgrna_dict = pickle.load(f)

base_editors = make_preset_base_editors()

formatted_exploded_sgrna_df = format_output(target_exons_df_with_sgrna_dict, base_editors)

formatted_exploded_sgrna_df.to_pickle("data/formatted_exploded_sgrna_df.pkl")