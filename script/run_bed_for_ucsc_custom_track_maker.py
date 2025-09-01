from  altex_be.bed_for_ucsc_custom_track_maker import format_sgrna_for_ucsc_custom_track
import pandas as pd
from pathlib import Path
import datetime

output_directory = Path("data/")
output_track_name = "altex_be_sgrnas"
exploded_sgrna_with_offtarget_info = pd.read_pickle("data/exploded_sgrna_df_with_offtarget_counts.pkl")

bed_df = format_sgrna_for_ucsc_custom_track(exploded_sgrna_with_offtarget_info)

output_path = output_directory / f"{output_track_name}_ucsc_custom_track.bed"
track_description: str = f"sgRNAs designed by altex-be on {datetime.datetime.now().strftime('%Y%m%d')}"

with open(output_path, "w") as f:
    track_header = f'track name="{output_track_name}" description="{track_description}" visibility=2 itemRgb="On"\n'
    f.write(track_header)
    bed_df.to_csv(f, sep="\t", header=False, index=False, lineterminator='\n')