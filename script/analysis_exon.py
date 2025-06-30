import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

pd.set_option('display.width', 200) 

data = pd.read_pickle('../data/exon_classification_with_additional_info.pkl')
print(data.head())
print(data.shape)

# - 全ての遺伝子について、バリアント数の分布を可視化

plt.figure(figsize=(8, 5))
variant_counts_per_gene = data.groupby("geneName")["variant_count"].first().values
sns.histplot(variant_counts_per_gene, bins=30, kde=False, color="skyblue")

plt.xlabel("Number of Transcript Variants per Gene")
plt.ylabel("Number of Genes")
plt.title("Distribution of Transcript Variant Counts per Gene")
plt.grid(True)
plt.yscale('log')
plt.tight_layout()
plt.show()


# - 遺伝子当たりの最大or最小エキソン数のヒストグラムを描画する
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


# - エキソン数がtop5の遺伝子の中身を確認する
top_5_exon_counts = data.nlargest(5, "max_exon_count")
print(top_5_exon_counts)


# - 最大エキソン数 vs バリアント数の散布図を作成するためのdfを作成
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=data,
    x="max_exon_count",
    y="variant_count",
    alpha=0.6,
    hue="coding",
    palette={"coding": "skyblue", "non-coding": "#D55E00"},
)

plt.xlabel("Max Exon Count per Gene")
plt.ylabel("Transcript Variant Count per Gene")
plt.title("Exon Count vs Transcript Diversity")
plt.grid(True)
plt.tight_layout()
plt.show()


# - skipped exon, unique exon, a3ss, a5ss, overlapを持つ遺伝子の数を計算
total_genes = data["geneName"].nunique()
print(f"全遺伝子数: {total_genes}")

# （skipped exon）を持つ遺伝子の数をカウント
genes_have_skipped_exon = data[data["exontype"].str.contains("skipped")]["geneName"].unique()
print(f"skipped exonを持つ遺伝子数: {len(genes_have_skipped_exon)}")
print(f"skipped exonを持たない遺伝子数: {total_genes - len(genes_have_skipped_exon)}")
print(f"skipped exonを持つ遺伝子の割合: {len(genes_have_skipped_exon) / total_genes:.2%}")

#  （unique exon）を持つ遺伝子の数をカウント
genes_have_unique_exon = data[data["exontype"].str.contains("unique")]["geneName"].unique()
print(f"unique exonを持つ遺伝子数: {len(genes_have_unique_exon)}")
print(f"unique exonを持たない遺伝子数: {total_genes - len(genes_have_unique_exon)}")
print(f"unique exonを持つ遺伝子の割合: {len(genes_have_unique_exon) / total_genes:.2%}")

# 少なくとも（skipped exon）と（unique exon）のいずれかを持つ遺伝子の数をカウント
genes_have_skipped_or_unique_exon = data[
    data["exontype"].str.contains("skipped|unique")]["geneName"].unique()
print(f"skipped exonまたはunique exonを持つ遺伝子数: {len(genes_have_skipped_or_unique_exon)}")
print(f"skipped exonとunique exonの両方を持たない遺伝子数: {total_genes - len(genes_have_skipped_or_unique_exon)}")
print(f"skipped exonまたはunique exonを持つ遺伝子の割合: {len(genes_have_skipped_or_unique_exon) / total_genes:.2%}")

# (a3ss alternative 3' splice site)を持つ遺伝子の数をカウント
genes_have_a3ss = data[data["exontype"].str.contains("a3ss")]["geneName"].unique()
print(f"a3ss alternative 3' splice siteを持つ遺伝子数: {len(genes_have_a3ss)}")
print(f"a3ss alternative 3' splice siteを持たない遺伝子数: {total_genes - len(genes_have_a3ss)}")
print(f"a3ss alternative 3' splice siteを持つ遺伝子の割合: {len(genes_have_a3ss) / total_genes:.2%}")

# (a5ss alternative 5' splice site)を持つ遺伝子の数をカウント
genes_have_a5ss = data[data["exontype"].str.contains("a5ss")]["geneName"].unique()
print(f"a5ss alternative 5' splice siteを持つ遺伝子数: {len(genes_have_a5ss)}")
print(f"a5ss alternative 5' splice siteを持たない遺伝子数: {total_genes - len(genes_have_a5ss)}")
print(f"a5ss alternative 5' splice siteを持つ遺伝子の割合: {len(genes_have_a5ss) / total_genes:.2%}")

# overlapping exonを持つ遺伝子の数をカウント
genes_have_overlapping_exon = data[data["exontype"].apply(lambda x: "overlapping" in x)]["geneName"].unique()
print(f"overlapping exonを持つ遺伝子数: {len(genes_have_overlapping_exon)}")
print(f"overlapping exonを持たない遺伝子数: {total_genes - len(genes_have_overlapping_exon)}")
print(f"overlapping exonを持つ遺伝子の割合: {len(genes_have_overlapping_exon) / total_genes:.2%}")

# split exonを持つ遺伝子の数をカウント
genes_have_intron_retention = data[data["exontype"].str.contains("intron_retention")]["geneName"].unique()
print(f"intron_retentionを持つ遺伝子数: {len(genes_have_intron_retention)}")
print(f"intron_retentionを持たない遺伝子数: {total_genes - len(genes_have_intron_retention)}")
print(f"intron_retentionを持つ遺伝子の割合: {len(genes_have_intron_retention) / total_genes:.2%}")

# outflame exonを持つ遺伝子の数をカウント
genes_have_outflame_exon = data[data["flame"].str.contains("out-flame")]["geneName"].unique()
print(f"outflame exonを持つ遺伝子数: {len(genes_have_outflame_exon)}")
print(f"outflame exonを持たない遺伝子数: {total_genes - len(genes_have_outflame_exon)}")
print(f"outflame exonを持つ遺伝子の割合: {len(genes_have_outflame_exon) / total_genes:.2%}")


# - 全遺伝子ではなく、skipped exonが生じうる可能性のある条件を満たすエキソンだけを分母としてみる

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
    data["exontype"].str.contains("skipped|unique")]["geneName"].unique()
print(f"skipped exonまたはunique exonを持つ遺伝子数: {len(genes_have_skipped_or_unique_exon)}")

# total genes数をカウント
total_genes = data["geneName"].nunique()
print(f"全遺伝子数: {total_genes}")

# skipped exonまたはunique exonが生じうる遺伝子を母数にしてSkipped exonを持つ遺伝子の割合を計算
print(f"skipped exonを持つ遺伝子の割合 (skipped exonまたはunique exonが生じうる遺伝子を母数とする): "
    f"{len(genes_have_skipped_or_unique_exon) / ((total_genes)-len(genes_with_one_exon_or_one_variant)):.2%}")

# 最大exonが一つまたはtranscript variantsが一つの遺伝子でskipped exonまたはunique exonを持つ遺伝子
genes_with_one_exon_or_one_variant_and_skipped_or_unique = data[
    (data["geneName"].isin(genes_with_one_exon_or_one_variant)) &
    (data["exontype"].str.contains("skipped|unique"))
]
print(f"最大エキソンが一つまたはtranscript variantsが一つの遺伝子でskipped exonまたはunique exonを持つ遺伝子数: "
    f"{len(genes_with_one_exon_or_one_variant_and_skipped_or_unique)}")