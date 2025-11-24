# %%
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from plotnine import*

# %%
pd.set_option('display.width', 200) 
pd.set_option('display.max_columns', 1000)
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
data2 = data[["geneName","exontype", "chrom", "exons", "coding","exonlengths", "frame", "exon_position", "cds_info"]]
exon_df = data2.explode(["exontype", "exons", "frame", "exonlengths", "exon_position", "cds_info"])
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
alternative_or_unique_alternative_exon = exon_df[exon_df["exontype"].isin(["alternative", "unique-alternative"])]
a5ss_long_exon = exon_df[exon_df["exontype"] == "a5ss-long"]
a3ss_long_exon = exon_df[exon_df["exontype"] == "a3ss-long"]
a5ss_short_exon = exon_df[exon_df["exontype"] == "a5ss-short"]
a3ss_short_exon = exon_df[exon_df["exontype"] == "a3ss-short"]
constitutive_exon = exon_df[exon_df["exontype"] == "constitutive"]


# %%
# cds_info のvalueごとに in-frame 割合を計算し棒グラフ化
def calc_in_frame_ratio(df):
    if len(df) == 0:
        return 0.0
    return (df["frame"] == "in-frame").sum() / len(df)

# 対象リスト
exon_sets = [
    ("alternative_exon", alternative_or_unique_alternative_exon),
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
        "cds_edge_exon_start_end": "utr_containing_exon",
        "cds_edge_exon_start": "utr_containing_exon",
        "cds_edge_exon_end": "utr_containing_exon",
        "utr_exon": "utr_containing_exon"
    })
    for cds_status in ["cds_exon", "utr_containing_exon"]:
        sub = df[df["cds_info_grouped"] == cds_status]
        total_count = sub.shape[0]
        in_frame_count = (sub["frame"] == "in-frame").sum()
        out_frame_count = (sub["frame"] == "out-frame").sum()
        ratio = calc_in_frame_ratio(sub)
        exon_count = sub.shape[0]
        results.append({
            "Category": name,
            "cording_status": cds_status,
            "total_exon_count": total_count,
            "in_frame_count": in_frame_count,
            "out_frame_count": out_frame_count,
            "InframeRatio": ratio * 100,
            "ExonCount": exon_count
        })

df_in_frame = pd.DataFrame(results)
print(df_in_frame)

# カテゴリの順序を定義
base_categories = ["constitutive_exon", "alternative_exon", "a5ss_long_exon", "a3ss_long_exon", "a5ss_short_exon", "a3ss_short_exon"]
full_categories = [f"{prefix}{cat}" for prefix in ["internal_", "first_", "last_"] for cat in base_categories]

df_in_frame["Category"] = pd.Categorical(df_in_frame["Category"], categories=["constitutive_exon", "alternative_exon", "a5ss_long_exon", "a3ss_long_exon", "a5ss_short_exon", "a3ss_short_exon"], ordered=True)
df_in_frame["cording_status"] = pd.Categorical(df_in_frame["cording_status"], categories=["cds_exon","utr_containing_exon"], ordered=True)

df_in_frame["label"] = df_in_frame.apply(
    lambda row: f"({int(row['ExonCount'])})", axis=1
)
exon_counts_df = df_in_frame[df_in_frame["cording_status"].isin(["cds_exon","utr_containing_exon"])]

for df in [exon_counts_df]:
    plot = (
        ggplot(df, aes(x="Category", y="InframeRatio", fill="cording_status")) +
        geom_bar(stat="identity", position="dodge") +
        labs(title="% of in-frame exon in splicing categories", x="Exon Category", y="% of in-frame exon") +
        scale_fill_manual(values={
            "cds_exon": "#E69F00", #orange
            "utr_containing_exon": "#009E73" #green
        }) +
        theme(
            axis_text_x=element_text(rotation=0, hjust=0.5, size=10),
            axis_title_y=element_text(size=10),
            axis_title_x=element_text(size=12),
            legend_title=element_text(size=10),
            legend_text=element_text(size=10),
            figure_size=(8,6)
        ) +
        scale_y_continuous(limits=(0, 100)) +
        coord_flip()
    )
display(plot)

plot.save("../data/exon_inframe_ratio_by_coding_status.png", dpi=600)


# %%
def summarize_coding_vs_noncoding_with_upstream(df: pd.DataFrame):
    """
    CDS alternative exon の out-frame のみ抽出し、
    coding / non-coding inclusion と upstream_alternative を同時に保持した pivot を作る
    """
    # explode して exon 単位に展開
    df_exp = df.explode(["exons", "exontype", "cds_info","frame", "structural_alternative"]).reset_index(drop=True)
    df_exp = df_exp[["geneName", "chrom", "exons", "coding", "cds_info", "frame", "exontype","structural_alternative"]]
    # exonsでgroupby して、cds_info を集約
    df_exp["cds_info_contains_cds_exon"] = df_exp.groupby(["geneName", "exons"])["cds_info"].transform(lambda x: x.str.contains("cds_exon").any())
    # CDS exon かつ alternative exon & out-frame のみ抽出
    df_alt_cds = df_exp[
        (df_exp["cds_info_contains_cds_exon"]) &
        (df_exp["exontype"].isin(["alternative", "unique-alternative"])) &
        (df_exp["frame"] == "out-frame")
    ]

    # coding / non-coding inclusion 集計
    coding_summary = (
        df_alt_cds.groupby(["geneName", "chrom", "exons", "coding", "structural_alternative"])
        .size()
        .reset_index(name="count")
    )
    print(coding_summary["coding"].value_counts())

    # pivot
    pivot = coding_summary.pivot_table(
        index=["geneName", "chrom","exons", "structural_alternative"],
        columns="coding",
        values="count",
        fill_value=0
    ).reset_index()

    pivot.columns.name = None

    # 列名統一
    cols = pivot.columns.tolist()
    if "coding" in cols:
        pivot = pivot.rename(columns={"coding": "coding_count"})
    if "non-coding" in cols:
        pivot = pivot.rename(columns={"non-coding": "noncoding_count"})

    # フラグ追加
    pivot["only_in_coding"] = (pivot.get("coding", pivot.get("coding_count", 0)) > 0) & \
                              (pivot.get("non-coding", pivot.get("noncoding_count", 0)) == 0)
    pivot["only_in_noncoding"] = (pivot.get("coding", pivot.get("coding_count", 0)) == 0) & \
                                 (pivot.get("non-coding", pivot.get("noncoding_count", 0)) > 0)

    return pivot
# 使用例
pivot = summarize_coding_vs_noncoding_with_upstream(data)

# %%
print(len(pivot))
print(len(pivot[(pivot["structural_alternative"] == True) & (pivot["only_in_coding"] == False) & (pivot["only_in_noncoding"] == False)]))
print(len(pivot[(pivot["structural_alternative"] == True) & (pivot["only_in_coding"] == True) & (pivot["only_in_noncoding"] == False)]))
print(len(pivot[(pivot["structural_alternative"] == True) & (pivot["only_in_noncoding"] == True) & (pivot["only_in_coding"] == False)]))
print(len(pivot[(pivot["structural_alternative"] == False) & (pivot["only_in_coding"] == True)]))
print(len(pivot[(pivot["structural_alternative"] == False) & (pivot["only_in_noncoding"] == True)]))
print(len(pivot[(pivot["structural_alternative"] == False) & (pivot["only_in_coding"] == False) & (pivot["only_in_noncoding"] == False)]))

print(pivot[(pivot["structural_alternative"] == False) & (pivot["only_in_coding"] == False) & (pivot["only_in_noncoding"] == False)].head(10))

# %%
plot_data = pd.DataFrame({
    "total_out_frame_exons": [len(pivot)],
    "structural_alternative_exons": [len(pivot[pivot["structural_alternative"] == True])],
    "only_in_coding_transcripts": [len(pivot[(pivot["only_in_coding"] == True) & (pivot["structural_alternative"] == False)])],
    "other_types": [len(pivot[(pivot["structural_alternative"] == False) & (pivot["only_in_coding"] == False)])]
})

# total を除いたデータを取得
plot_data_no_total = plot_data.drop(columns=["total_out_frame_exons"])

# 各値の割合を計算
sizes = plot_data_no_total.iloc[0]
percentages = sizes / sizes.sum() * 100

# 値が大きい順に並べ替え
sorted_data = pd.DataFrame({
    "Category": plot_data_no_total.columns,
    "Size": sizes,
    "Percentage": percentages
}).sort_values(by="Size", ascending=False)

# 円グラフを作成
labels = sorted_data["Category"]  # 並べ替えたラベル
sizes = sorted_data["Size"]       # 並べ替えた値

# パーセンテージと値を表示する関数
def autopct_with_count(pct, all_vals):
    absolute = int(round(pct / 100. * sum(all_vals)))
    return f"{pct:.1f}%\n(n={absolute})"

# プロット
# 新しいカラーパレット（オレンジ、淡いオレンジ、グレー）
colors = ['#E69F00', '#F0E68C', '#D3D3D3']  # オレンジ、淡いオレンジ、グレー

# プロット
plt.figure(figsize=(8, 8))
plt.pie(
    sizes, 
    labels=labels, 
    autopct=lambda pct: autopct_with_count(pct, sizes), 
    startangle=90,  # 時計回りに開始
    counterclock=False,  # 時計回りに描画
    colors=colors,  # 新しいカラーパレット
    textprops={'fontsize': 15},
)
plt.title("Distribution of Out-Frame Exons \n CDS-Alternative-exon (Percentage of Total)", fontsize=16, pad=20)
plt.axis('equal')  # 円を丸く表示
plt.savefig("../data/out_frame_exon_distribution_sorted.png", dpi=600)
plt.show()

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
colors = ['#009E73', '#66CDAA', '#D3D3D3']
plt.figure(figsize=(6, 6))
plt.pie(
    sizes, 
    labels=labels, 
    colors=colors, 
    autopct=lambda pct: autopct_with_count(pct, sizes), 
    startangle=90,
    textprops={'fontsize': '14'},
    wedgeprops={'edgecolor': 'black', 'linewidth': 1}  # 枠線の設定
)
plt.title('Proportion of genes that can theoretically have splicing events (mm39)', pad=20)
plt.axis('equal')  # 円を丸く表示
plt.savefig('../data/mm39/proportion_of_genes_mm39_with_borders.png', dpi=300)
plt.show()


# %%
# 関数を呼び出して結果を取得
result_df_based_on_theoretically_valid_conditions = count_genes_by_exon_type(genes_not_one_exon_or_one_variant, exon_conditions)
# 結果を表示
print(result_df_based_on_theoretically_valid_conditions)

result_df_based_on_theoretically_valid_conditions = result_df_based_on_theoretically_valid_conditions[
    result_df_based_on_theoretically_valid_conditions["Exon Type"].isin([
        "Alternative or unique-alternative exon",
        "Alternative 3' splice site-long",
        "Alternative 5' splice site-long",
        "Alternative 3' splice site-short",
        "Alternative 5' splice site-short",
        "Overlapping exon",
        "Intron retention"
    ])
]

print(result_df_based_on_theoretically_valid_conditions)

result_df_based_on_theoretically_valid_conditions["Exon Type"] = result_df_based_on_theoretically_valid_conditions["Exon Type"].replace({
    "Overlapping exon": "Overlapping exon",
    "Intron retention": "Intron retention",
    "Alternative 3' splice site-short": "A3SS-short",
    "Alternative 5' splice site-short": "A5SS-short",
    "Alternative 3' splice site-long": "A3SS-long",
    "Alternative 5' splice site-long": "A5SS-long",
    "Alternative or unique-alternative exon": "Alternative exon"
})

# Exon Type列の順序を明示的に指定
result_df_based_on_theoretically_valid_conditions["Exon Type"] = pd.Categorical(
    result_df_based_on_theoretically_valid_conditions["Exon Type"],
    categories=[
        "Overlapping exon", 
        "Intron retention",
        "A5SS-short",
        "A3SS-short",
        "A5SS-long",
        "A3SS-long",
        "Alternative exon"
    ],  # 希望する順序を指定
    ordered=True
)
result_df_based_on_theoretically_valid_conditions = result_df_based_on_theoretically_valid_conditions.dropna(subset=["Exon Type"])
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
            "A3SS-long": "#E69F00",  # オレンジ
            "A5SS-long": "#E69F00",  # オレンジ
            "A3SS-short": "#56B4E9",  # 青
            "A5SS-short": "#56B4E9",  # 青
            "Overlapping exon": "#56B4E9",  # 青
            "Intron retention": "#56B4E9",  # 青
            "Genes have targetable splicing events": "#009E73"  # 緑
        }
    ) +
    labs(
        title="% of Genes containing splicing events \n(Denominator: Genes have multiple isoforms)",
        x="Exon Type",
        y="% of Genes have each splicing event",
        fill="Category"
    ) +
    theme(
        axis_text_x=element_text(size=11, rotation=0, hjust=0.5),
        axis_text_y=element_text(size=11, rotation=0, hjust=1),
        figure_size=(8, 6),
    ) + 
    scale_y_continuous(limits=(0, 100)) +
    coord_flip()
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


