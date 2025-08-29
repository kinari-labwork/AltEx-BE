import pandas as pd
from pathlib import Path

def format_sgrna_for_ucsc_custom_track(
    sgrna_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Purpose : 最終出力のsgRNA情報をUCSCカスタムトラック用のBED形式に変換する
    Parameters:
        sgrna_df (pd.DataFrame): offtarget までの情報を含むsgRNA情報のDataFrame。
    Return : pd.DataFrame 12 bedに修正された DataFrame
    """
    # score 列に100を超える値が入ることがあるため、100でクリップする
    sgrna_df["pam+20bp_exact_match"] = sgrna_df["pam+20bp_exact_match_count"].apply(lambda x: 100 if x > 100 else x)

    bed_df = pd.DataFrame()
    bed_df["chrom"] = sgrna_df["chrom"]
    bed_df["chromStart"] = sgrna_df["sgrna_start_in_genome"]
    bed_df["chromEnd"] = sgrna_df["sgrna_end_in_genome"]
    bed_df["name"] = sgrna_df["geneName"] + "_" + sgrna_df["site_type"] + "_" + sgrna_df["base_editor_name"]  + "_" + sgrna_df["uuid"].str[:8]
    bed_df["score"] = sgrna_df["pam+20bp_exact_match_count"]
    bed_df["strand"] = sgrna_df["strand"]
    bed_df["thickStart"] = bed_df["chromStart"] # 別に必要ないが、9bed にするために追加
    bed_df["thickEnd"] = bed_df["chromEnd"]

    # Add RGB color based on base editor type
    color_map = {"abe": "255,0,0", "cbe": "0,0,255"} # ABE: red, CBE: blue
    bed_df["itemRgb"] = sgrna_df["base_editor_type"].map(color_map)
    
    # Reorder columns for BED9 format
    bed_df = bed_df[["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"]]
    return bed_df