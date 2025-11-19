# %%
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from plotnine import*

# %%
pd.set_option('display.width', 200) 
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 1000)

# %%
data = pd.read_pickle('../data/mm39/classified_refflat_all_genes.pkl')
print(data.head())
print(data.shape)

# %% [markdown]
# - 全ての遺伝子について、バリアント数の分布を可視化

# %% [markdown]
# - 遺伝子当たりの最大or最小エキソン数のヒストグラムを描画する

# %%
def max_min_exon_count_annotator(data: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームに、最大および最小のエクソン数を追加する。
    Parameters:
        data: pd.DataFrame, refFlatのデータフレーム
    Returns:
        pd.DataFrame, 最大および最小のエクソン数を追加したrefFlatのデータフレーム
    """
    maximum_exon_counts = data.groupby("geneName")["exonCount"].max().reset_index()
    maximum_exon_counts.columns = ["geneName", "max_exon_count"]
    minimum_exon_counts = data.groupby("geneName")["exonCount"].min().reset_index()
    minimum_exon_counts.columns = ["geneName", "min_exon_count"]
    data = data.merge(maximum_exon_counts, on="geneName", how="left")
    data = data.merge(minimum_exon_counts, on="geneName", how="left")
    return data

# %%
# 最大および最小のエクソン数を追加
data = max_min_exon_count_annotator(data)

# %%
def add_variant_count(data: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose:
        refFlatのデータフレームに、variant_count列を追加する。
    Parameters:
        data: pd.DataFrame, refFlatのデータフレーム
    Returns:
        pd.DataFrame, variant_count列を追加したrefFlatのデータフレーム
    """
    # geneNameごとにvariant_countを計算
    variant_counts = data.groupby("geneName").size().reset_index(name="variant_count")
    # データフレームにvariant_count列を追加
    data = data.merge(variant_counts, on="geneName", how="left")
    return data

# 関数を使用してvariant_count列を追加
data = add_variant_count(data)

# %%
print(data.columns)

# %%
maxcount = data.groupby("geneName")["max_exon_count"].first().values
mincount = data.groupby("geneName")["min_exon_count"].first().values
ange_exon_counts = maxcount - mincount

plt.figure(figsize= (8,5))
sns.histplot(maxcount, bins=50, kde=False, color="skyblue")
plt.xlabel("Number of maximum exons per Gene")
plt.ylabel("Number of Genes")
plt.title("Distribution of Maximum Exon Counts per Gene")
plt.grid(True)
plt.yscale('log')
plt.tight_layout()
plt.show()

plt.figure(figsize= (8,5))
sns.histplot(data["min_exon_count"], bins=50, kde=False, color="#D55E00")
plt.xlabel("Number of minimum exons per Gene")
plt.ylabel("Number of Genes")
plt.title("Distribution of Minimum Exon Counts per Gene")
plt.grid(True)
plt.yscale('log')
plt.tight_layout()
plt.show()

range_exon_counts = data.groupby("geneName")["exonCount"].apply(lambda x: max(x) - min(x)).values

plt.figure(figsize=(8, 5))
sns.histplot(range_exon_counts, bins=50, kde=False, color="lightgreen")
plt.xlabel("Range of Exon Counts per Gene")
plt.ylabel("Number of Genes")
plt.title("Distribution of Exon Count Range per Gene")
plt.grid(True)
plt.yscale('log')
plt.tight_layout()
plt.show()

# %% [markdown]
# - エキソン数がtop5の遺伝子の中身を確認する

# %%
top_5_exon_counts = data.nlargest(5, "max_exon_count")
print(top_5_exon_counts)

# %% [markdown]
# - 最大エキソン数 vs バリアント数の散布図を作成するためのdfを作成

# %% [markdown]
# - skipped exon, unique exon, a3ss, a5ss, overlapを持つ遺伝子の数を計算

# %%
total_genes = data["geneName"].nunique()
print(f"全遺伝子数: {total_genes}")

# （skipped exon）を持つ遺伝子の数をカウント
genes_have_alternative_exon = data[data["exontype"].apply(lambda x: "alternative" in x)]["geneName"].unique()

print(f"alternative exonを持つ遺伝子数: {len(genes_have_alternative_exon)}")
print(f"skipped exonを持たない遺伝子数: {total_genes - len(genes_have_alternative_exon)}")
print(f"skipped exonを持つ遺伝子の割合: {len(genes_have_alternative_exon) / total_genes:.2%}")

#  （unique exon）を持つ遺伝子の数をカウント
genes_have_unique_exon = data[data["exontype"].apply(lambda x: "unique-alternative" in x)]["geneName"].unique()
print(f"unique-alternative exonを持つ遺伝子数: {len(genes_have_unique_exon)}")
print(f"unique-alternative exonを持たない遺伝子数: {total_genes - len(genes_have_unique_exon)}")
print(f"unique-alternative exonを持つ遺伝子の割合: {len(genes_have_unique_exon) / total_genes:.2%}")

# 少なくとも（skipped exon）と（unique exon）のいずれかを持つ遺伝子の数をカウント
genes_have_skipped_or_unique_exon = data[
    data["exontype"].apply(lambda x: "alternative" in x or "unique-alternative" in x)]["geneName"].unique()
print(f"skipped exonまたはunique exonを持つ遺伝子数: {len(genes_have_skipped_or_unique_exon)}")
print(f"skipped exonとunique exonの両方を持たない遺伝子数: {total_genes - len(genes_have_skipped_or_unique_exon)}")
print(f"skipped exonまたはunique exonを持つ遺伝子の割合: {len(genes_have_skipped_or_unique_exon) / total_genes:.2%}")

# (a3ss alternative 3' splice site)を持つ遺伝子の数をカウント
genes_have_a3ss = data[data["exontype"].apply(lambda x: ("a3ss-short" in x or "a3ss-long" in x))]["geneName"].unique()
print(f"a3ss alternative 3' splice siteを持つ遺伝子数: {len(genes_have_a3ss)}")
print(f"a3ss alternative 3' splice siteを持たない遺伝子数: {total_genes - len(genes_have_a3ss)}")
print(f"a3ss alternative 3' splice siteを持つ遺伝子の割合: {len(genes_have_a3ss) / total_genes:.2%}")

# (a5ss alternative 5' splice site)を持つ遺伝子の数をカウント
genes_have_a5ss = data[data["exontype"].apply(lambda x: ("a5ss-short" in x or "a5ss-long" in x))]["geneName"].unique()
print(f"a5ss alternative 5' splice siteを持つ遺伝子数: {len(genes_have_a5ss)}")
print(f"a5ss alternative 5' splice siteを持たない遺伝子数: {total_genes - len(genes_have_a5ss)}")
print(f"a5ss alternative 5' splice siteを持つ遺伝子の割合: {len(genes_have_a5ss) / total_genes:.2%}")

# overlapping exonを持つ遺伝子の数をカウント
genes_have_overlapping_exon = data[data["exontype"].apply(lambda x: "overlap" in x)]["geneName"].unique()
print(f"overlapping exonを持つ遺伝子数: {len(genes_have_overlapping_exon)}")
print(f"overlapping exonを持たない遺伝子数: {total_genes - len(genes_have_overlapping_exon)}")
print(f"overlapping exonを持つ遺伝子の割合: {len(genes_have_overlapping_exon) / total_genes:.2%}")

# split exonを持つ遺伝子の数をカウント
genes_have_intron_retention = data[data["exontype"].apply(lambda x: "intron_retention" in x)]["geneName"].unique()
print(f"intron_retentionを持つ遺伝子数: {len(genes_have_intron_retention)}")
print(f"intron_retentionを持たない遺伝子数: {total_genes - len(genes_have_intron_retention)}")
print(f"intron_retentionを持つ遺伝子の割合: {len(genes_have_intron_retention) / total_genes:.2%}")

# %%
print(data.head(2))

# %%
data2 = data[["geneName","exontype", "chrom", "exons", "coding","exonlengths", "flame", "exon_position", "cds_info"]]
exon_df = data2.explode(["exontype", "exons", "flame", "exonlengths", "exon_position", "cds_info"])
exon_df = exon_df.drop_duplicates(subset=["chrom","exons"])
print(exon_df.head())

# %%
print(exon_df["coding"].value_counts())

# %%
print(exon_df["cds_info"].value_counts())

# %%
target_exon = exon_df[exon_df["exontype"].isin(["alternative", "unique-alternative", "a5ss-long", "a3ss-long", "a5ss-short", "a3ss-short"])]
internal_target_exon = target_exon[target_exon["exon_position"] == "internal"]
first_target_exon = target_exon[target_exon["exon_position"] == "first"]
last_target_exon = target_exon[target_exon["exon_position"] == "last"]

internal_unique_exon = internal_target_exon[internal_target_exon["exontype"] == "unique-alternative"]
first_unique_exon = first_target_exon[first_target_exon["exontype"] == "unique-alternative"]
last_unique_exon = last_target_exon[last_target_exon["exontype"] == "unique-alternative"]

internal_skipped_exon = internal_target_exon[internal_target_exon["exontype"] == "alternative"]
first_skipped_exon = first_target_exon[first_target_exon["exontype"] == "alternative"]
last_skipped_exon = last_target_exon[last_target_exon["exontype"] == "alternative"]

internal_a5ss_long_exon = internal_target_exon[internal_target_exon["exontype"] == "a5ss-long"]
first_a5ss_long_exon = first_target_exon[first_target_exon["exontype"] == "a5ss-long"]
last_a5ss_long_exon = last_target_exon[last_target_exon["exontype"] == "a5ss-long"]

internal_a3ss_long_exon = internal_target_exon[internal_target_exon["exontype"] == "a3ss-long"]
first_a3ss_long_exon = first_target_exon[first_target_exon["exontype"] == "a3ss-long"]
last_a3ss_long_exon = last_target_exon[last_target_exon["exontype"] == "a3ss-long"] 

internal_a5ss_short_exon = internal_target_exon[internal_target_exon["exontype"] == "a5ss-short"]
first_a5ss_short_exon = first_target_exon[first_target_exon["exontype"] == "a5ss-short"]
last_a5ss_short_exon = last_target_exon[last_target_exon["exontype"] == "a5ss-short"]

internal_a3ss_short_exon = internal_target_exon[internal_target_exon["exontype"] == "a3ss-short"]
first_a3ss_short_exon = first_target_exon[first_target_exon["exontype"] == "a3ss-short"]
last_a3ss_short_exon = last_target_exon[last_target_exon["exontype"] == "a3ss-short"]

print(internal_skipped_exon["cds_info"].value_counts())

# %%
alternative_exon = exon_df[exon_df["exontype"] == "alternative"]
unique_alternative_exon = exon_df[exon_df["exontype"] == "unique-alternative"]
a5ss_long_exon = exon_df[exon_df["exontype"] == "a5ss-long"]
a3ss_long_exon = exon_df[exon_df["exontype"] == "a3ss-long"]
a5ss_short_exon = exon_df[exon_df["exontype"] == "a5ss-short"]
a3ss_short_exon = exon_df[exon_df["exontype"] == "a3ss-short"]
constitutive_exon = exon_df[exon_df["exontype"] == "constitutive"]


# %%
# cds_info のvalueごとに in-flame 割合を計算し棒グラフ化
def calc_in_flame_ratio(df):
    if len(df) == 0:
        return 0.0
    return (df["flame"] == "in-flame").sum() / len(df)

# 対象リスト
exon_sets = [
    ("unique_alternative_exon", unique_alternative_exon),
    ("alternative_exon", alternative_exon),
    ("a5ss_long_exon", a5ss_long_exon),
    ("a3ss_long_exon", a3ss_long_exon),
    ("a5ss_short_exon", a5ss_short_exon),
    ("a3ss_short_exon", a3ss_short_exon),
    ("constitutive_exon", constitutive_exon)
]

# 結果格納用
results = []
for name, df in exon_sets:
    # "cds_edge_exon_start_end", "cds_edge_exon_start", "cds_edge_exon_end" をまとめて "cds_edge_exon" に
    df = df.copy()
    df["cds_info_grouped"] = df["cds_info"].replace({
        "cds_edge_exon_start_end": "cds_edge_exon",
        "cds_edge_exon_start": "cds_edge_exon",
        "cds_edge_exon_end": "cds_edge_exon"
    })
    for cds_status in ["cds_exon", "utr_exon", "cds_edge_exon"]:
        sub = df[df["cds_info_grouped"] == cds_status]
        total_count = sub.shape[0]
        in_flame_count = (sub["flame"] == "in-flame").sum()
        out_flame_count = (sub["flame"] == "out-flame").sum()
        ratio = calc_in_flame_ratio(sub)
        exon_count = sub.shape[0]
        results.append({
            "Category": name,
            "cording_status": cds_status,
            "total_exon_count": total_count,
            "in_flame_count": in_flame_count,
            "out_flame_count": out_flame_count,
            "InFlameRatio": ratio * 100,
            "ExonCount": exon_count
        })

df_in_flame = pd.DataFrame(results)
print(df_in_flame.head())

# カテゴリの順序を定義
base_categories = ["constitutive_exon", "unique_alternative_exon", "alternative_exon", "a5ss_long_exon", "a3ss_long_exon", "a5ss_short_exon", "a3ss_short_exon"]
full_categories = [f"{prefix}{cat}" for prefix in ["internal_", "first_", "last_"] for cat in base_categories]

df_in_flame["Category"] = pd.Categorical(df_in_flame["Category"], categories=["constitutive_exon", "unique_alternative_exon", "alternative_exon", "a5ss_long_exon", "a3ss_long_exon", "a5ss_short_exon", "a3ss_short_exon"], ordered=True)
df_in_flame["cording_status"] = pd.Categorical(df_in_flame["cording_status"], categories=["cds_exon","cds_edge_exon","utr_exon"], ordered=True)

df_in_flame["label"] = df_in_flame.apply(
    lambda row: f"({int(row['ExonCount'])})", axis=1
)
exon_counts_df = df_in_flame[df_in_flame["cording_status"].isin(["cds_exon","utr_exon","cds_edge_exon"])]

for df in [exon_counts_df]:
    plot = (
        ggplot(df, aes(x="Category", y="InFlameRatio", fill="cording_status")) +
        geom_bar(stat="identity", position="dodge") +
        geom_text(
            aes(label="label"),
            position=position_dodge(width=0.9),
            va="bottom",
            size=11,
            color="black"
        ) +
        labs(title="Percentage of in-flame exon in splicing categories", x="Exon Category", y="percentage of in-flame exon (%)") +
        coord_cartesian(ylim=(0, 100)) +
        scale_fill_manual(values={
            "cds_exon": "#56B4E9", #blue
            "utr_exon": "#009E73", #green
            "cds_edge_exon": "#E69F00", #orange
        }) +
        theme(
            axis_text_x=element_text(rotation=90, hjust=0.5, size=15),
            axis_title_y=element_text(size=15),
            legend_title=element_text(size=15),
            legend_text=element_text(size=12),
            figure_size=(17,10)
        )
    )
    display(plot)

# %%
# coding列が "coding" または "non-coding" の場合の in-flame割合を計算し棒グラフ化

def calc_in_flame_ratio(df):
    if len(df) == 0:
        return 0.0
    return (df["flame"] == "in-flame").sum() / len(df)

# 対象リスト
exon_sets = [
    ("internal_unique_alternative_exon", internal_unique_exon),
    ("first_unique_alternative_exon", first_unique_exon),
    ("last_unique_alternative_exon", last_unique_exon),
    ("internal_alternative_exon", internal_skipped_exon),
    ("first_alternative_exon", first_skipped_exon),
    ("last_alternative_exon", last_skipped_exon),
    ("internal_a5ss_long_exon", internal_a5ss_long_exon),
    ("first_a5ss_long_exon", first_a5ss_long_exon),
    ("last_a5ss_long_exon", last_a5ss_long_exon),
    ("internal_a3ss_long_exon", internal_a3ss_long_exon),
    ("first_a3ss_long_exon", first_a3ss_long_exon),
    ("last_a3ss_long_exon", last_a3ss_long_exon),
]

# 結果格納用
results = []
for name, df in exon_sets:
    for coding_status in ["coding", "non-coding"]:
        sub = df[df["coding"] == coding_status]
        ratio = calc_in_flame_ratio(sub)
        results.append({"Category": name, "Coding": coding_status, "InFlameRatio": ratio * 100})

df_in_flame = pd.DataFrame(results)

# 棒グラフ
plot = (
    ggplot(df_in_flame, aes(x="Category", y="InFlameRatio", fill="Coding")) +
    geom_bar(stat="identity", position="dodge") +
    geom_text(
        aes(label="InFlameRatio.round(1).astype(str) + '%'"),
        position=position_dodge(width=0.9),
        va="bottom",
        size=9,
        color="black"
    ) +
    labs(title="Percentage of in-flame exon in splicing categories", x="Exon Category", y="percentage of in-flame exon (%)") +
    coord_cartesian(ylim=(0, 100)) +
    theme(axis_text_x=element_text(rotation=45, hjust=1), figure_size=(12,5))
)
plot.save("../data/in_flame_ratio.png", dpi=600)
display(plot)

# %%
# 各項目のin-flame割合をまとめる
in_flame_ratios = {
    "internal_unique_exon": calc_in_flame_ratio(internal_unique_exon),
    "first_unique_exon": calc_in_flame_ratio(first_unique_exon),
    "last_unique_exon": calc_in_flame_ratio(last_unique_exon),
    "internal_skipped_exon": calc_in_flame_ratio(internal_skipped_exon),
    "first_skipped_exon": calc_in_flame_ratio(first_skipped_exon),
    "last_skipped_exon": calc_in_flame_ratio(last_skipped_exon),
    "internal_a5ss_long_exon": calc_in_flame_ratio(internal_a5ss_long_exon),
    "first_a5ss_long_exon": calc_in_flame_ratio(first_a5ss_long_exon),
    "last_a5ss_long_exon": calc_in_flame_ratio(last_a5ss_long_exon),
    "internal_a3ss_long_exon": calc_in_flame_ratio(internal_a3ss_long_exon),
    "first_a3ss_long_exon": calc_in_flame_ratio(first_a3ss_long_exon),
    "last_a3ss_long_exon": calc_in_flame_ratio(last_a3ss_long_exon),
}

# データフレーム化
df_in_flame = pd.DataFrame({
    "Category": list(in_flame_ratios.keys()),
    "InFlameRatio": [v * 100 for v in in_flame_ratios.values()]
})

# 棒グラフ
plot = (
    ggplot(df_in_flame, aes(x="Category", y="InFlameRatio")) +
    geom_bar(stat="identity", fill="#56B4E9") +
    geom_text(aes(label="InFlameRatio.round(2).astype(str) + '%'"), va="bottom", size=9, color="black") +
    labs(title="Percentage of in-flame Exons", x="Exon Category", y="Percentage of in-flame exons (%)") +
    coord_cartesian(ylim=(0, 100)) +
    theme(axis_text_x=element_text(rotation=45, hjust=1), figure_size=(10,5))
)
display(plot)

# %%
def count_genes_by_exon_type(data, exon_conditions):
    """
    各エキソンタイプごとに遺伝子数をカウントし、結果を出力する関数。

    Parameters:
        data (pd.DataFrame): エキソンデータを含むデータフレーム。
        exon_conditions (list of tuple): 各エキソンタイプのラベルと条件を指定するリスト。
            例: [("skipped exon", lambda x: "skipped" in x), ...]

    Returns:
        pd.DataFrame: 各エキソンタイプごとの遺伝子数、割合をまとめたデータフレーム。
    """
    total_genes = data["geneName"].nunique()
    results = []

    for label, condition in exon_conditions:
        genes = data[data["exontype"].apply(condition)]["geneName"].unique()
        count = len(genes)
        not_count = total_genes - count
        ratio = (count / total_genes)*100
        results.append({
            "Exon Type": label,
            "Genes with Exon": count,
            "Genes without Exon": not_count,
            "Percentage": f"{ratio:.2f}"
        })

    return pd.DataFrame(results)

# 各エキソンタイプの条件を定義
exon_conditions = [
    ("Alternative exon", lambda x: "alternative" in x),
    ("Unique-Alternative exon", lambda x: "unique-alternative" in x),
    ("Alternative or unique-alternative exon", lambda x: "alternative" in x or "unique-alternative" in x),
    ("Alternative 3' splice site-long", lambda x: "a3ss-long" in x),
    ("Alternative 5' splice site-long", lambda x: "a5ss-long" in x),
    ("Alternative 3' splice site-short", lambda x: "a3ss-short" in x),
    ("Alternative 5' splice site-short", lambda x: "a5ss-short" in x),
    ("Overlapping exon", lambda x: "overlap" in x),
    ("Intron retention", lambda x: "intron_retention" in x),
    ("Genes have targetable splicing events", lambda x: "alternative" in x or "unique-alternative" in x or "a3ss-long" in x  or "a5ss-long" in x)
]

# 関数を呼び出して結果を取得
result_df = count_genes_by_exon_type(data, exon_conditions)

# %%
# 結果を表示
print(result_df)

# %%
# Exon Type列の順序を明示的に指定
result_df["Exon Type"] = pd.Categorical(
    result_df["Exon Type"],
    categories=[
        "Alternative exon", 
        "Unique-Alternative exon", 
        "Alternative or unique-alternative exon", 
        "Alternative 3' splice site-long", 
        "Alternative 5' splice site-long", 
        "Alternative 3' splice site-short",
        "Alternative 5' splice site-short",
        "Overlapping exon", 
        "Intron retention",
        "Genes have targetable splicing events"
    ],  # 希望する順序を指定
    ordered=True
)

# データをlong形式に変換
plot_df = pd.melt(
    result_df,
    id_vars=["Exon Type", "Percentage"],
    value_vars=["Genes with Exon", "Genes without Exon"],
    var_name="Category",
    value_name="Count"
)
# カテゴリ名をわかりやすく
plot_df["Category"] = plot_df["Category"].replace({
    "Genes with Exon": "Genes With Splicing Event",
    "Genes without Exon": "Genes Without Splicing Event"
})

plot_df["Category"] = pd.Categorical(
    plot_df["Category"],
    categories=["Genes Without Splicing Event", "Genes With Splicing Event"],  # 希望する順序を指定
    ordered=True
)

# Percentage列を数値型に変換
plot_df["Percentage"] = plot_df["Percentage"].str.rstrip('%').astype(float)
print(plot_df)



# プロット
plot = (
    ggplot(plot_df, aes(x="Exon Type", y="Count", fill="Category")) +
    geom_bar(stat="identity", position="stack") +
    geom_text(
        aes(label="Percentage"),
        data=plot_df[plot_df["Category"] == "Genes With Splicing Event"],  # "With Exon" のみラベルを表示
        position="stack",
        va="bottom",  # バーの上に配置
        size=10,
        color="black",
        format_string="{:.1f}%"  # パーセンテージを表示
    ) +
    scale_fill_manual(
        values={
            "Genes With Splicing Event": "#56B4E9",  # 青
            "Genes Without Splicing Event": "#E69F00"  # オレンジ
        }
    ) +
    labs(
        title="Percentage of Genes by Exon Type (Denominator: Total Genes)",
        x="Exon Type",
        y="Number of Genes",
        fill="Category"
    ) +
    theme(
        axis_text_x=element_text(size=11, rotation=90, hjust=1),
        figure_size=(10, 6),
    )
)

display(plot)  # Jupyter Notebookでの表示用

# %% [markdown]
# - 全遺伝子ではなく、splicing eventが生じうる可能性のある条件を満たすエキソンだけを分母としてみる

# %%
genes_not_one_exon_or_one_variant = data[
    (data["max_exon_count"] != 1) & (data["variant_count"] != 1)
    ]
genes_one_exon = data[data["max_exon_count"] == 1]
genes_not_one_exon_but_one_variant = data[(data["variant_count"] == 1) & (data["max_exon_count"] != 1)]

def autopct_with_count(pct, all_vals):
    absolute = int(round(pct / 100. * sum(all_vals)))
    return f"{pct:.1f}%\n(n={absolute})"

# 原理的にsplice event が生じうる遺伝子の割合を円グラフで表示
total_genes = data["geneName"].nunique()
valid_genes = genes_not_one_exon_or_one_variant["geneName"].nunique()
genes_one_exon_count = genes_one_exon["geneName"].nunique()
genes_not_one_exon_but_one_variant_count = genes_not_one_exon_but_one_variant["geneName"].nunique()

labels = ['Genes have multiple isoforms', 'Genes have a single isoform', 'Genes have one exon']
sizes = [valid_genes, genes_not_one_exon_but_one_variant_count, genes_one_exon_count,]
colors = ['#E69F00', '#009E73', '#56B4E9']
plt.figure(figsize=(6, 6))
plt.pie(
    sizes, 
    labels=labels, 
    colors=colors, 
    autopct=lambda pct: autopct_with_count(pct, sizes), 
    startangle=90,
    textprops={'fontsize': '14'}
    )
plt.title('Proportion of genes that can theoretically have splicing events (mm39)', pad=20)
plt.savefig('../data/mm39/proportion_of_genes_mm39.png', dpi=300)
plt.axis('equal')  # 円を丸く表示
plt.show()


# %%
# 関数を呼び出して結果を取得
result_df_based_on_theoretically_valid_conditions = count_genes_by_exon_type(genes_not_one_exon_or_one_variant, exon_conditions)
# 結果を表示
# Exon Type列の順序を明示的に指定
result_df_based_on_theoretically_valid_conditions["Exon Type"] = pd.Categorical(
    result_df_based_on_theoretically_valid_conditions["Exon Type"],
    categories=[
        "Alternative exon", 
        "Unique-Alternative exon", 
        "Alternative or unique-alternative exon", 
        "Alternative 3' splice site-long", 
        "Alternative 5' splice site-long", 
        "Alternative 3' splice site-short",
        "Alternative 5' splice site-short",
        "Overlapping exon", 
        "Intron retention",
        "Genes have targetable splicing events"
    ],  # 希望する順序を指定
    ordered=True
)

# データをlong形式に変換
plot_df = pd.melt(
    result_df_based_on_theoretically_valid_conditions,
    id_vars=["Exon Type", "Percentage"],
    value_vars=["Genes with Exon", "Genes without Exon"],
    var_name="Category",
    value_name="Count"
)
# カテゴリ名をわかりやすく
plot_df["Category"] = plot_df["Category"].replace({
    "Genes with Exon": "Genes With Splicing Event",
    "Genes without Exon": "Genes Without Splicing Event"
})

plot_df["Category"] = pd.Categorical(
    plot_df["Category"],
    categories=["Genes Without Splicing Event", "Genes With Splicing Event"],  # 希望する順序を指定
    ordered=True
)

# Percentage列を数値型に変換
plot_df["Percentage"] = plot_df["Percentage"].str.rstrip('%').astype(float)

plotdf = plot_df[plot_df["Category"] == "Genes With Splicing Event"]
# プロット
plot = (
    ggplot(plot_df, aes(x="Exon Type", y="Percentage", fill="Exon Type")) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(
        values={
            "Alternative exon": "#E69F00",  # オレンジ
            "Unique-Alternative exon": "#E69F00",  # オレンジ
            "Alternative or unique-alternative exon": "#E69F00",  # オレンジ
            "Alternative 3' splice site-long": "#E69F00",  # オレンジ
            "Alternative 5' splice site-long": "#E69F00",  # オレンジ
            "Alternative 3' splice site-short": "#56B4E9",  # 青
            "Alternative 5' splice site-short": "#56B4E9",  # 青
            "Overlapping exon": "#56B4E9",  # 青
            "Intron retention": "#56B4E9",  # 青
            "Genes have targetable splicing events": "#009E73"  # 緑
        }
    ) +
    labs(
        title="Percentage of Genes containing splicing events \n(Denominator: Genes have multiple isoforms)",
        x="Exon Type",
        y="Percentage of Genes \nhave each splicing event (%)",
        fill="Category"
    ) +
    theme(
        axis_text_x=element_text(size=11, rotation=90, hjust=0.5),
        figure_size=(8, 6),
    )
)

display(plot)  # Jupyter Notebookでの表示用
plot.save("../data/mm39/percentage_genes_has_splicing_events_mm39.png", dpi=300)

# %%
# 最大エキソンが1つの遺伝子をカウント
genes_with_one_exon = data[data["max_exon_count"] == 1]["geneName"].unique()
print(f"最大エキソンが1つの遺伝子数: {len(genes_with_one_exon)}")
# transcript variantsが1つの遺伝子をカウント
genes_with_one_variant = data[data["variant_count"] == 1]["geneName"].unique()  
print(f"transcript variantsが1つの遺伝子数: {len(genes_with_one_variant)}")



# %%
# エキソンが一つまたはtranscript variantsが一つの遺伝子をカウント
genes_with_one_exon_or_one_variant = data[
    (data["max_exon_count"] == 1) | (data["variant_count"] == 1)
]["geneName"].unique()
print(f"最大エキソンが一つまたはtranscript variantsが一つの遺伝子数: {len(genes_with_one_exon_or_one_variant)}")

# skipped exon or unique exonをもつ遺伝子をカウント
genes_have_skipped_or_unique_exon = data[
    data["exontype"].apply(lambda x: "skipped" in x or "unique" in x)]["geneName"].unique()

# total genes数をカウント
total_genes = data["geneName"].nunique()
print(f"全遺伝子数: {total_genes}")

# skipped exonまたはunique exonが生じうる遺伝子を母数にしてSkipped exonを持つ遺伝子の割合を計算
print(f"skipped exonを持つ遺伝子の割合 (skipped exonまたはunique exonが生じうる遺伝子を母数とする): "
    f"{len(genes_have_skipped_or_unique_exon) / ((total_genes)-len(genes_with_one_exon_or_one_variant)):.2%}")

# 最大exonが一つまたはtranscript variantsが一つの遺伝子でskipped exonまたはunique exonを持つ遺伝子
genes_with_one_exon_or_one_variant_and_skipped_or_unique = data[
    (data["geneName"].isin(genes_with_one_exon_or_one_variant)) &
    (data["exontype"].apply(lambda x: "skipped" in x or "unique" in x))
]
print(f"最大エキソンが一つまたはtranscript variantsが一つの遺伝子でskipped exonまたはunique exonを持つ遺伝子数: "
    f"{len(genes_with_one_exon_or_one_variant_and_skipped_or_unique)}")


