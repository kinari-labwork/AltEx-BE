# %%
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from plotnine import*

# %%
pd.set_option('display.width', 200) 

# %%
data = pd.read_pickle('../data/classified_exon_refflat.pkl')
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
genes_have_skipped_exon = data[data["exontype"].apply(lambda x: "skipped" in x)]["geneName"].unique()

print(f"skipped exonを持つ遺伝子数: {len(genes_have_skipped_exon)}")
print(f"skipped exonを持たない遺伝子数: {total_genes - len(genes_have_skipped_exon)}")
print(f"skipped exonを持つ遺伝子の割合: {len(genes_have_skipped_exon) / total_genes:.2%}")

#  （unique exon）を持つ遺伝子の数をカウント
genes_have_unique_exon = data[data["exontype"].apply(lambda x: "unique" in x)]["geneName"].unique()
print(f"unique exonを持つ遺伝子数: {len(genes_have_unique_exon)}")
print(f"unique exonを持たない遺伝子数: {total_genes - len(genes_have_unique_exon)}")
print(f"unique exonを持つ遺伝子の割合: {len(genes_have_unique_exon) / total_genes:.2%}")

# 少なくとも（skipped exon）と（unique exon）のいずれかを持つ遺伝子の数をカウント
genes_have_skipped_or_unique_exon = data[
    data["exontype"].apply(lambda x: "skipped" in x or "unique" in x)]["geneName"].unique()
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
    ("skipped exon", lambda x: "skipped" in x),
    ("unique exon", lambda x: "unique" in x),
    ("skipped or unique exon", lambda x: "skipped" in x or "unique" in x),
    ("alternative 3' splice site-long", lambda x: "a3ss-long" in x),
    ("alternative 5' splice site-long", lambda x: "a5ss-long" in x),
    ("alternative 3' splice site-short", lambda x: "a3ss-short" in x),
    ("alternative 5' splice site-short", lambda x: "a5ss-short" in x),
    ("overlapping exon", lambda x: "overlap" in x),
    ("intron retention", lambda x: "intron_retention" in x),
    ("targetable gene", lambda x: "skipped" in x or "unique" in x or "a3ss-long" in x  or "a5ss-long" in x)
]

# 関数を呼び出して結果を取得
result_df = count_genes_by_exon_type(data, exon_conditions)

# 結果を表示
print(result_df)

# %%
# Exon Type列の順序を明示的に指定
result_df["Exon Type"] = pd.Categorical(
    result_df["Exon Type"],
    categories=[
        "skipped exon", 
        "unique exon", 
        "skipped or unique exon", 
        "alternative 3' splice site-long", 
        "alternative 5' splice site-long", 
        "alternative 3' splice site-short",
        "alternative 5' splice site-short",
        "overlapping exon", 
        "intron retention",
        "targetable gene"
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
    "Genes with Exon": "Genes With Applicable Type Exon",
    "Genes without Exon": "Genes Without Applicable Type Exon"
})

plot_df["Category"] = pd.Categorical(
    plot_df["Category"],
    categories=["Genes Without Applicable Type Exon", "Genes With Applicable Type Exon"],  # 希望する順序を指定
    ordered=True
)

# Percentage列を数値型に変換
plot_df["Percentage"] = plot_df["Percentage"].str.rstrip('%').astype(float)

# プロット
plot = (
    ggplot(plot_df, aes(x="Exon Type", y="Count", fill="Category")) +
    geom_bar(stat="identity", position="stack") +
    geom_text(
        aes(label="Percentage"),
        data=plot_df[plot_df["Category"] == "Genes With Applicable Type Exon"],  # "With Exon" のみラベルを表示
        position="stack",
        va="bottom",  # バーの上に配置
        size=10,
        color="black",
        format_string="{:.1f}%"  # パーセンテージを表示
    ) +
    scale_fill_manual(
        values={
            "Genes With Applicable Type Exon": "#56B4E9",  # 青
            "Genes Without Applicable Type Exon": "#E69F00"  # オレンジ
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
# - 全遺伝子ではなく、skipped exonが生じうる可能性のある条件を満たすエキソンだけを分母としてみる

# %%
genes_not_one_exon_or_one_variant = data[
    (data["max_exon_count"] != 1) & (data["variant_count"] != 1)
]

# 関数を呼び出して結果を取得
result_df_based_on_theoretically_valid_conditions = count_genes_by_exon_type(genes_not_one_exon_or_one_variant, exon_conditions)
# 結果を表示
print(result_df_based_on_theoretically_valid_conditions)

# %%
# Exon Type列の順序を明示的に指定
result_df_based_on_theoretically_valid_conditions["Exon Type"] = pd.Categorical(
    result_df_based_on_theoretically_valid_conditions["Exon Type"],
    categories=[
        "skipped exon", 
        "unique exon", 
        "skipped or unique exon", 
        "alternative 3' splice site-long", 
        "alternative 5' splice site-long", 
        "alternative 3' splice site-short",
        "alternative 5' splice site-short",
        "overlapping exon", 
        "intron retention",
        "targetable gene"
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
    "Genes with Exon": "Genes With Applicable Type Exon",
    "Genes without Exon": "Genes Without Applicable Type Exon"
})

plot_df["Category"] = pd.Categorical(
    plot_df["Category"],
    categories=["Genes Without Applicable Type Exon", "Genes With Applicable Type Exon"],  # 希望する順序を指定
    ordered=True
)

# Percentage列を数値型に変換
plot_df["Percentage"] = plot_df["Percentage"].str.rstrip('%').astype(float)

# プロット
plot = (
    ggplot(plot_df, aes(x="Exon Type", y="Count", fill="Category")) +
    geom_bar(stat="identity", position="stack") +
    geom_text(
        aes(label="Percentage"),
        data=plot_df[plot_df["Category"] == "Genes With Applicable Type Exon"],  # "With Exon" のみラベルを表示
        position="stack",
        va="bottom",  # バーの上に配置
        size=10,
        color="black",
        format_string="{:.1f}%"  # パーセンテージを表示
    ) +
    scale_fill_manual(
        values={
            "Genes With Applicable Type Exon": "#56B4E9",  # 青
            "Genes Without Applicable Type Exon": "#E69F00"  # オレンジ
        }
    ) +
    labs(
        title="Percentage of Genes by Exon Type (Denominator: theoretically Valid Genes)",
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

# %%
# 最大エキソンが1つの遺伝子をカウント
genes_with_one_exon = data[data["max_exon_count"] == 1]["geneName"].unique()
print(f"最大エキソンが1つの遺伝子数: {len(genes_with_one_exon)}")
# transcript variantsが1つの遺伝子をカウント
genes_with_one_variant = data[data["variant_count"] == 1]["geneName"].unique()  
print(f"transcript variantsが1つの遺伝子数: {len(genes_with_one_variant)}")

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


