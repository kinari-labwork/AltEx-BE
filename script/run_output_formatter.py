import pandas as pd
from pathlib import Path
import pickle
from altex_aid.output_formatter import format_output
from altex_aid.sgrna_designer import BaseEditor, make_preset_base_editors

with open("target_exon_with_sgrna_dict.pkl", "rb") as f:
    target_exons_df_with_sgrna_dict = pickle.load(f)

exploded_classified_refflat = pd.read_csv("exploded_classified_refflat.csv")

base_editors = make_preset_base_editors()

formatted_exploded_sgrna_df = format_output(target_exons_df_with_sgrna_dict, exploded_classified_refflat, base_editors)
