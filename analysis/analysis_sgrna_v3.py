# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %%

pd.set_option('display.width', 100000000)
pd.set_option("display.max_columns", None)  # すべての列を表示
pd.set_option("display.max_colwidth", None)
pd.set_option('display.max_rows', 5000)  # 表示する行数の上限を設定
P = print

# %%
res_df = pd.read_csv("../data/mm39/altexbe_20260326_164828/202603261736_mm39_sgrnas_designed_by_altex-be_table.csv")

# %%

def calc_offtarget_stat(ref_df, base_editors, output_csv_path="../data/mm39/offtarget_statistics.csv"):
    """
    Calculate comprehensive off-target statistics for each base editor.
    
    Parameters:
        ref_df: DataFrame with sgRNA data
        base_editors: List of base editor names
        output_csv_path: Path to save the statistics CSV
        
    Returns:
        DataFrame with statistics for all base editors
    """
    stats_list = []
    
    for base_editor in base_editors:
        be_df = ref_df[ref_df["base_editor_name"] == base_editor]
        
        # 20bp statistics
        off_target_20mer = be_df["pam+20bp_exact_match_count"]
        # 12bp statistics
        off_target_12mer = be_df["pam+12bp_exact_match_count"]
        
        stats_dict = {
            "base_editor": base_editor,
            "designed_gRNA_count": len(off_target_20mer),
            # 20bp
            "20bp_mean": off_target_20mer.mean(),
            "20bp_std": off_target_20mer.std(),
            "20bp_median": off_target_20mer.median(),
            "20bp_min": off_target_20mer.min(),
            "20bp_max": off_target_20mer.max(),
            "20bp_Q25": off_target_20mer.quantile(0.25),
            "20bp_Q50": off_target_20mer.quantile(0.50),
            "20bp_Q75": off_target_20mer.quantile(0.75),
            "20bp_Q95": off_target_20mer.quantile(0.95),
            "20bp_vector": off_target_20mer.tolist(),
            # 12bp
            "12bp_mean": off_target_12mer.mean(),
            "12bp_std": off_target_12mer.std(),
            "12bp_median": off_target_12mer.median(),
            "12bp_min": off_target_12mer.min(),
            "12bp_max": off_target_12mer.max(),
            "12bp_Q25": off_target_12mer.quantile(0.25),
            "12bp_Q50": off_target_12mer.quantile(0.50),
            "12bp_Q75": off_target_12mer.quantile(0.75),
            "12bp_Q95": off_target_12mer.quantile(0.95),
            "12bp_vector": off_target_12mer.tolist(),
        }
        stats_list.append(stats_dict)
    
    stats_df = pd.DataFrame(stats_list)
    stat_df_to_save = stats_df.drop(columns=["20bp_vector", "12bp_vector"])
    stat_df_to_save.to_csv(output_csv_path, index=False)
    print(f"✓ Statistics saved to: {output_csv_path}")
    
    return stats_df

base_editors = ["target_aid_ngg", "be4max_ngg", "abe8e_ngg", "target_aid_ng", "be4max_ng", "abe8e_ng", "SpRY_BE4max"]

# Calculate and export off-target statistics
stats_df = calc_offtarget_stat(res_df, base_editors)
print("\n" + "="*80)
print(stats_df.to_string())
print("="*80)

# %%
# Add PAM type and base editor type for visualization
def add_metadata(stats_df):
    """Add PAM type and editor type for better visualization."""
    def get_pam_type(editor_name):
        if "ngg" in editor_name.lower():
            return "NGG"
        elif "ng" in editor_name.lower():
            return "NG"
        elif "nnn" in editor_name.lower():
            return "NNN"
        return "Other"
    
    def get_editor_type(editor_name):
        if "abe" in editor_name.lower():
            return "ABE"
        else:
            return "CBE"
    
    stats_df = stats_df.copy()
    stats_df["pam_type"] = stats_df["base_editor"].apply(get_pam_type)
    stats_df["editor_type"] = stats_df["base_editor"].apply(get_editor_type)
    return stats_df

stats_df = add_metadata(stats_df)

# %%
# Plot: Box plot for 20bp and 12bp off-target counts by base editor

plot_20bp_df = stats_df[["base_editor", "20bp_vector"]].explode("20bp_vector")
plot_20bp_df["20bp_vector"] = pd.to_numeric(plot_20bp_df["20bp_vector"], errors="coerce")
plot_20bp_df = plot_20bp_df.dropna(subset=["20bp_vector"])

plot_12bp_df = stats_df[["base_editor", "12bp_vector"]].explode("12bp_vector")
plot_12bp_df["12bp_vector"] = pd.to_numeric(plot_12bp_df["12bp_vector"], errors="coerce")
plot_12bp_df = plot_12bp_df.dropna(subset=["12bp_vector"])

plt.figure(figsize=(12, 6))

# --- 20bp ---
plt.subplot(1, 2, 1)
sns.violinplot(
    x="base_editor",
    y="20bp_vector",
    data=plot_20bp_df,
    inner=None,  
    cut=0,              # ← 外側に変な裾を伸ばさない
    scale="width"       # ← 比較しやすい
)
plt.yscale("log")
plt.title("Distribution of 20bp Off-target Counts by Base Editor")
plt.xticks(rotation=30, ha="right")

# --- 12bp ---
plt.subplot(1, 2, 2)
sns.violinplot(
    x="base_editor",
    y="12bp_vector",
    data=plot_12bp_df,
    inner=None,
    cut=0,
    scale="width"
)
plt.yscale("log")
plt.title("Distribution of 12bp Off-target Counts by Base Editor")
plt.xticks(rotation=30, ha="right")

plt.tight_layout()

# %%
# Summary table: Compare by PAM type
print("\n" + "="*80)
print("Summary by PAM Type:")
print("="*80)
summary_by_pam = stats_df.groupby("pam_type")[["designed_gRNA_count", "20bp_mean", "12bp_mean", "12bp_median"]].agg(["mean", "min", "max"])
print(summary_by_pam.round(2))
