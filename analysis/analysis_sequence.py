# %% [markdown]
# 取得した配列の解析

# %%
import pandas as pd
from plotnine import *
import logomaker
import matplotlib.pyplot as plt


# %%
pd.set_option('display.width', 100000)
pd.set_option("display.max_colwidth", 10000)
P = print

# %%
data = pd.read_pickle("../data/target_exons_with_acceptor_and_donor_sequence.pkl")

# %%
print(data[data["strand"] == "+"].head(10))

# %%
internal_plus_exons = (data["exon_position"] == "internal") & (data["strand"] == "+")
internal_minus_exons = (data["exon_position"] == "internal") & (data["strand"] == "-")
internal_exons = internal_plus_exons | internal_minus_exons
first_plus_exons = (data["exon_position"] == "first") & (data["strand"] == "+")
first_minus_exons = (data["exon_position"] == "first") & (data["strand"] == "-")
first_exons = first_plus_exons | first_minus_exons
last_plus_exons = (data["exon_position"] == "last") & (data["strand"] == "+")
last_minus_exons = (data["exon_position"] == "last") & (data["strand"] == "-")
last_exons = last_plus_exons | last_minus_exons


print("\n--- general information ---\n")
print("Number of exons in the data:", len(data))
print("Number of internal exons:", internal_exons.sum())
print("Number of first exons:", first_exons.sum())
print("Number of last exons:", last_exons.sum())

# %%
target_acceptor_regions = data["exon_position"].isin(["internal", "last"])
target_donor_regions = data["exon_position"].isin(["internal", "first"])
is_AG = data.loc[target_acceptor_regions, "acceptor_sequence"].str[23:25] == "AG"
is_ag = data.loc[target_acceptor_regions, "acceptor_sequence"].str[23:25] == "ag"
is_GT = data.loc[target_donor_regions, "donor_sequence"].str[25:27] == "GT"
is_gt = data.loc[target_donor_regions, "donor_sequence"].str[25:27] == "gt"

print("\n--- integrated information ---\n")
print("Number of target exons with acceptor sequence:", target_acceptor_regions.sum())
print("Number of target exons with donor sequence:", target_donor_regions.sum())
print("Number of target exons with AG acceptor sequence:", is_AG.sum())
print("Number of target exons with ag acceptor sequence:", is_ag.sum())
print("Number of target exons with GT donor sequence:", is_GT.sum())
print("Number of target exons with gt donor sequence:", is_gt.sum())
print("Number of target exons with AG or ag acceptor sequence:", (is_AG | is_ag).sum())
print("Number of target exons with GT or gt donor sequence:", (is_GT | is_gt).sum())
print("Percentage of target exons with AG or ag acceptor sequence:", (is_AG | is_ag).sum() / target_acceptor_regions.sum() * 100)
print("Percentage of target exons with GT or gt donor sequence:", (is_GT | is_gt).sum() / target_donor_regions.sum() * 100)

# %%
is_AG_plus = data.loc[internal_plus_exons, "acceptor_sequence"].str[23:25] == "AG"
is_GT_plus = data.loc[internal_plus_exons, "donor_sequence"].str[25:27] == "GT"
is_ag_plus = data.loc[internal_plus_exons, "acceptor_sequence"].str[23:25] == "ag"
is_gt_plus = data.loc[internal_plus_exons, "donor_sequence"].str[25:27] == "gt"

print("\n--- internal plus exons ---\n")
print("Number of internal plus exons:", internal_plus_exons.sum())
print("Number of internal plus exons with AG acceptor sequence:", is_AG_plus.sum())
print("Number of internal plus exons with GT donor sequence:", is_GT_plus.sum())
print("Number of internal plus exons with ag acceptor sequence:", is_ag_plus.sum())
print("Number of internal plus exons with gt donor sequence:", is_gt_plus.sum())
print("Number of internal plus exons with AG or ag acceptor sequence:", (is_AG_plus | is_ag_plus).sum())
print("Number of internal plus exons with GT or gt donor sequence:", (is_GT_plus | is_gt_plus).sum())
print("Percentage of internal plus exons with AG acceptor sequence:", (is_AG_plus| is_ag_plus).sum()/ internal_plus_exons.sum() * 100)
print("Percentage of internal plus exons with GT donor sequence:", (is_GT_plus | is_gt_plus).sum() / internal_plus_exons.sum() * 100)

# %%
is_AG_minus = data.loc[internal_minus_exons, "acceptor_sequence"].str[23:25] == "AG"
is_GT_minus = data.loc[internal_minus_exons, "donor_sequence"].str[25:27] == "GT"
is_ag_minus = data.loc[internal_minus_exons, "acceptor_sequence"].str[23:25] == "ag"
is_gt_minus = data.loc[internal_minus_exons, "donor_sequence"].str[25:27] == "gt"

print("\n--- Internal Minus Exons ---\n")
print("Number of internal minus exons:", internal_minus_exons.sum())
print("Number of internal minus exons with AG acceptor sequence:", is_AG_minus.sum())
print("Number of internal minus exons with GT donor sequence:", is_GT_minus.sum())
print("Number of internal minus exons with ag acceptor sequence:", is_ag_minus.sum())
print("Number of internal minus exons with gt donor sequence:", is_gt_minus.sum())
print("Number of internal minus exons with AG or ag acceptor sequence:", (is_AG_minus | is_ag_minus).sum())
print("Number of internal minus exons with GT or gt donor sequence:", (is_GT_minus | is_gt_minus).sum())
print("Percentage of internal minus exons with AG acceptor sequence:", (is_AG_minus | is_ag_minus).sum() / internal_minus_exons.sum() * 100)
print("Percentage of internal minus exons with GT donor sequence:", (is_GT_minus | is_gt_minus).sum() / internal_minus_exons.sum() * 100)

# %%
is_AG_plus = data.loc[first_plus_exons, "acceptor_sequence"].str[23:25] == "AG"
is_GT_plus = data.loc[first_plus_exons, "donor_sequence"].str[25:27] == "GT"
is_ag_plus = data.loc[first_plus_exons, "acceptor_sequence"].str[23:25] == "ag"
is_gt_plus = data.loc[first_plus_exons, "donor_sequence"].str[25:27] == "gt"

print("\n--- first plus exons ---\n")
print("Number of first plus exons:", first_plus_exons.sum())
print("Number of first plus exons with AG acceptor sequence:", is_AG_plus.sum())
print("Number of first plus exons with GT donor sequence:", is_GT_plus.sum())
print("Number of first plus exons with ag acceptor sequence:", is_ag_plus.sum())
print("Number of first plus exons with gt donor sequence:", is_gt_plus.sum())
print("Number of first plus exons with AG or ag acceptor sequence:", (is_AG_plus | is_ag_plus).sum())
print("Number of first plus exons with GT or gt donor sequence:", (is_GT_plus | is_gt_plus).sum())
print("Percentage of first plus exons with AG acceptor sequence:", (is_AG_plus | is_ag_plus).sum() / first_plus_exons.sum() * 100)
print("Percentage of first plus exons with GT donor sequence:", (is_GT_plus | is_gt_plus).sum() / first_plus_exons.sum() * 100) 

# %%
is_AG_minus = data.loc[first_minus_exons, "acceptor_sequence"].str[23:25] == "AG"
is_GT_minus = data.loc[first_minus_exons, "donor_sequence"].str[25:27] == "GT"
is_ag_minus = data.loc[first_minus_exons, "acceptor_sequence"].str[23:25] == "ag"
is_gt_minus = data.loc[first_minus_exons, "donor_sequence"].str[25:27] == "gt"

print("\n--- first minus exons ---\n")
print("Number of first minus exons:", first_minus_exons.sum())
print("Number of first minus exons with AG acceptor sequence:", is_AG_minus.sum())
print("Number of first minus exons with GT donor sequence:", is_GT_minus.sum())
print("Number of first minus exons with ag acceptor sequence:", is_ag_minus.sum())
print("Number of first minus exons with gt donor sequence:", is_gt_minus.sum())
print("Number of first minus exons with AG or ag acceptor sequence:", (is_AG_minus | is_ag_minus).sum())
print("Number of first minus exons with GT or gt donor sequence:", (is_GT_minus | is_gt_minus).sum())
print("Percentage of first minus exons with AG acceptor sequence:", (is_AG_minus | is_ag_minus).sum() / first_minus_exons.sum() * 100)
print("Percentage of first minus exons with GT donor sequence:", (is_GT_minus | is_gt_minus).sum() / first_minus_exons.sum() * 100)

# %%
is_AG_plus = data.loc[last_plus_exons, "acceptor_sequence"].str[23:25] == "AG"
is_GT_plus = data.loc[last_plus_exons, "donor_sequence"].str[25:27] == "GT"
is_ag_plus = data.loc[last_plus_exons, "acceptor_sequence"].str[23:25] == "ag"
is_gt_plus = data.loc[last_plus_exons, "donor_sequence"].str[25:27] == "gt"
print("\n--- last plus exons ---\n")
print("Number of last plus exons:", last_plus_exons.sum())
print("Number of last plus exons with AG acceptor sequence:", is_AG_plus.sum())
print("Number of last plus exons with GT donor sequence:", is_GT_plus.sum())
print("Number of last plus exons with ag acceptor sequence:", is_ag_plus.sum())
print("Number of last plus exons with gt donor sequence:", is_gt_plus.sum())
print("Number of last plus exons with AG or ag acceptor sequence:", (is_AG_plus | is_ag_plus).sum())
print("Number of last plus exons with GT or gt donor sequence:", (is_GT_plus | is_gt_plus).sum())
print("Percentage of last plus exons with AG acceptor sequence:", (is_AG_plus | is_ag_plus).sum() / last_plus_exons.sum() * 100)
print("Percentage of last plus exons with GT donor sequence:", (is_GT_plus | is_gt_plus).sum() / last_plus_exons.sum() * 100)


# %%
is_AG_minus = data.loc[last_minus_exons, "acceptor_sequence"].str[23:25] == "AG"
is_GT_minus = data.loc[last_minus_exons, "donor_sequence"].str[25:27] == "GT"
is_ag_minus = data.loc[last_minus_exons, "acceptor_sequence"].str[23:25] == "ag"
is_gt_minus = data.loc[last_minus_exons, "donor_sequence"].str[25:27] == "gt"
print("\n--- last minus exons ---\n")
print("Number of last minus exons:", last_minus_exons.sum())
print("Number of last minus exons with AG acceptor sequence:", is_AG_minus.sum())
print("Number of last minus exons with GT donor sequence:", is_GT_minus.sum())
print("Number of last minus exons with ag acceptor sequence:", is_ag_minus.sum())
print("Number of last minus exons with gt donor sequence:", is_gt_minus.sum())
print("Number of last minus exons with AG or ag acceptor sequence:", (is_AG_minus | is_ag_minus).sum())
print("Number of last minus exons with GT or gt donor sequence:", (is_GT_minus | is_gt_minus).sum())
print("Percentage of last minus exons with AG acceptor sequence:", (is_AG_minus | is_ag_minus).sum() / last_minus_exons.sum() * 100)
print("Percentage of last minus exons with GT donor sequence:", (is_GT_minus | is_gt_minus).sum() / last_minus_exons.sum() * 100)

# %%
# 配列データを取得
internal_exons = data.loc[internal_exons,"acceptor_sequence"].dropna().str.upper()  # 欠損値を除去し、大文字に変換
first_exons = data.loc[first_exons, "acceptor_sequence"].dropna().str.upper()  # 欠損値を除去し、大文字に変換
last_exons = data.loc[last_exons, "acceptor_sequence"].dropna().str.upper()  # 欠損値を除去し、大文字に変換

seq_list = [internal_exons, first_exons, last_exons]

for i, seq in enumerate(seq_list):
    frequency_matrix = logomaker.alignment_to_matrix(seq.tolist(), to_type="probability")
    plt.figure(figsize=(15, 5))
    logo = logomaker.Logo(frequency_matrix)
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.set_ylabel("Frequency")
    logo.ax.set_xlabel("Position in Sequence")
    custom_ticks = list(range(0, frequency_matrix.shape[0], 10)) + [23, 24]  # 10刻み + 23, 24
    custom_ticks = sorted(set(custom_ticks))  # 重複を除去してソート
    logo.ax.set_xticks(custom_ticks)
    logo.ax.set_xticklabels(custom_ticks, rotation=45)
    plt.title(f"Splice acceptor sequence logo of {['Internal', 'First', 'Last'][i]} Exons")

    # シーケンスロゴを表示
    plt.show()

# %%
# 配列データを取得
internal_exons = (data["exon_position"] == "internal") & (data["strand"].isin(["+", "-"]))
internal_exons = data.loc[internal_exons,"donor_sequence"].dropna().str.upper()  # 欠損値を除去し、大文字に変換
first_exons = (data["exon_position"] == "first") & (data["strand"].isin(["+", "-"]))
first_exons = data.loc[first_exons, "donor_sequence"].dropna().str.upper()  # 欠損値を除去し、大文字に変換
last_exons = (data["exon_position"] == "last") & (data["strand"].isin(["+", "-"]))
last_exons = data.loc[last_exons, "donor_sequence"].dropna().str.upper()  # 欠損値を除去し、大文字に変換

seq_list = [internal_exons, first_exons, last_exons]

for i, seq in enumerate(seq_list):
    frequency_matrix = logomaker.alignment_to_matrix(seq.tolist(), to_type="probability")
    plt.figure(figsize=(15, 5))
    logo = logomaker.Logo(frequency_matrix)
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.set_ylabel("Frequency")
    logo.ax.set_xlabel("Position in Sequence")
    custom_ticks = list(range(0, frequency_matrix.shape[0], 10)) + [25, 26]  # 10刻み + 25, 26
    custom_ticks = sorted(set(custom_ticks))  # 重複を除去してソート
    logo.ax.set_xticks(custom_ticks)
    logo.ax.set_xticklabels(custom_ticks, rotation=45)
    plt.title(f"Splice donor sequence logo of {['Internal', 'First', 'Last'][i]} Exons")

    # シーケンスロゴを表示
    plt.show()


