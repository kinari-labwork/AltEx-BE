import pandas as pd
from pathlib import Path
import logging
import datetime
from . import logging_config # noqa: F401

def format_sgrna_for_ucsc_custom_track(
    sgrna_df: pd.DataFrame,
    output_directory: Path,
    track_name: str = "sgRNAs",
) -> None:
    """
    Converts the final sgRNA DataFrame to a BED file for UCSC custom track.

    Args:
        sgrna_df (pd.DataFrame): The final DataFrame with sgRNA information.
        output_directory (Path): The directory to save the BED file.
        track_name (str, optional): The name of the track. Defaults to "sgRNAs".
        track_description (str, optional): The description of the track. Defaults to "sgRNAs designed by Altex-BE".
    """

    bed_df = pd.DataFrame()
    bed_df["chrom"] = sgrna_df["chrom"]
    bed_df["chromStart"] = sgrna_df["sgrna_start_in_genome"]
    bed_df["chromEnd"] = sgrna_df["sgrna_end_in_genome"]
    bed_df["name"] = sgrna_df["geneName"] + "_" + sgrna_df["exontype"] + "_" + sgrna_df["base_editor_name"]  + "_" + sgrna_df["uuid"].str[:8]
    bed_df["score"] = sgrna_df["pam+20bp_exact_match"]
    bed_df["strand"] = sgrna_df["strand"]
    bed_df["thickStart"] = bed_df["sgrna_start_in_genome"]
    bed_df["thickEnd"] = bed_df["sgrna_end_in_genome"]

    # Add RGB color based on base editor type
    color_map = {"abe": "255,0,0", "cbe": "0,0,255"} # ABE: red, CBE: blue
    bed_df["itemRgb"] = sgrna_df["base_editor_type"].map(color_map)

    # Reorder columns for BED9 format
    bed_df = bed_df[["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"]]

    output_path = output_directory / f"{track_name}_ucsc_custom_track.bed"
    track_description: str = f"sgRNAs designed by Altex-BE on {datetime.datetime.now().strftime('%Y%m%d')}"

    with open(output_path, "w") as f:
        track_header = f'track name="{track_name}" description="{track_description}" visibility=2 itemRgb="On"\n'
        f.write(track_header)
        bed_df.to_csv(f, sep="\t", header=False, index=False, lineterminator='\n')

    logging.info(f"UCSC custom track file saved to: {output_path}")

