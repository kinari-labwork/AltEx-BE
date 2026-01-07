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

data = max_min_exon_count_annotator(data)

# %%
def annotate_variant_counts(data: pd.DataFrame) -> pd.DataFrame:
    """ 各遺伝子について、transcript variant の数をカウントする関数。 """ 
    variant_counts = data.groupby("geneName").size().reset_index(name="variant_count") 
    data = data.merge(variant_counts, on="geneName", how="left")
    return data

data = annotate_variant_counts(data)

# %%
data_filtered = data[["geneName","exontype", "chrom", "name","cdsStart", "cdsEnd","exons", "coding","exonlengths", "frame", "exon_position", "cds_info", "is_outside_common_exon_space", "max_exon_count", "min_exon_count", "variant_count"]]
exon_df = data_filtered.explode(["exontype", "exons", "frame", "exonlengths", "exon_position", "cds_info", "is_outside_common_exon_space"])
exon_df = exon_df.drop_duplicates(subset=["chrom","exons"])

# %%
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
base_categories = ["alternative_exon", "a5ss_long_exon", "a3ss_long_exon", "a5ss_short_exon", "a3ss_short_exon"]
full_categories = [f"{prefix}{cat}" for prefix in ["internal_", "first_", "last_"] for cat in base_categories]

df_in_frame["Category"] = pd.Categorical(df_in_frame["Category"], categories=["alternative_exon", "a5ss_long_exon", "a3ss_long_exon", "a5ss_short_exon", "a3ss_short_exon"], ordered=True)
df_in_frame = df_in_frame.dropna(subset=["Category"])
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
        scale_x_discrete(
            limits=[
                "alternative_exon",
                "a5ss_short_exon",
                "a3ss_short_exon",
                "a5ss_long_exon",
                "a3ss_long_exon",
            ],
            labels={
            "alternative_exon": "Alternative exon",
            "a5ss_long_exon": "A5SS-long",
            "a3ss_long_exon": "A3SS-long",
            "a5ss_short_exon": "A5SS-short",
            "a3ss_short_exon": "A3SS-short",
        }) +
        scale_y_continuous(limits=(0, 100)) +
        coord_flip()
    )
display(plot)

plot.save("../data/exon_inframe_ratio_by_coding_status.png", dpi=600)
# %%
# alternative かつ 少なくとも 1回 CDS exon となる exon を再分類し、 in-frame 割合を計算する
def extract_alternative_coding_exons(df_exp: pd.DataFrame) -> pd.DataFrame:
    """ alternative かつ 少なくとも 1回 CDS exon となるexonを抽出する関数。 """
    df_exp["cds_exon_any_coding"] = df_exp.groupby(["geneName", "exons"])["cds_info"].transform(
    lambda x: (x == "cds_exon").any()
    )
    df_exp = df_exp[
    (df_exp["exontype"].isin(["alternative", "unique-alternative"])) &
    (df_exp["cds_exon_any_coding"])
    ]
    return df_exp
# alternative かつ、CDS exon がなぜin-frame割合が低いのかを調査する
def annotate_variant_counts(data: pd.DataFrame) -> pd.DataFrame: 
    """ 各遺伝子について、coding variant と non-coding variant の数をカウントする関数。 """ 
    coding_counts = data[data["name"].str.startswith("NM")].groupby("geneName").size().reset_index(name="coding_variant_count") 
    noncoding_counts = data[data["name"].str.startswith("NR")].groupby("geneName").size().reset_index(name="noncoding_variant_count") 
    merged_counts = coding_counts.merge(noncoding_counts, on="geneName", how="outer").fillna(0) 
    data = data.merge(merged_counts, on="geneName", how="right") 
    return data


def annotate_exon_skipping(df_exp: pd.DataFrame) -> pd.DataFrame:
    """
    exon ごとに
    - skipped_only_in_coding
    - skipped_only_in_noncoding
    - skipped_both_coding_and_noncoding
    を判定
    """

    # 初期化
    df_exp["skipped_in_coding"] = False
    df_exp["skipped_only_in_noncoding"] = False

    # そこまで時間はかからないが、明らかに非効率であるので、将来的に改善を行うべき
    for gene, df_g in df_exp.groupby("geneName"):
        coding_total = df_g["coding_variant_count"].iloc[0]
        noncoding_total = df_g["noncoding_variant_count"].iloc[0]

        for exon, df_e in df_g.groupby("exons"):
            # transcript 単位で inclusion を数える
            coding_included = (
                df_e[df_e["name"].str.startswith("NM")]["name"].nunique()
            )
            noncoding_included = (
                df_e[df_e["name"].str.startswith("NR")]["name"].nunique()
            )

            skipped_in_coding = coding_included < coding_total
            skipped_in_noncoding = noncoding_included < noncoding_total

            if skipped_in_coding:
                flag = "skipped_in_coding"
            elif not skipped_in_coding and skipped_in_noncoding:
                flag = "skipped_only_in_noncoding"
            else:
                raise ValueError("Unexpected condition")

            df_exp.loc[
                (df_exp["geneName"] == gene) &
                (df_exp["exons"] == exon),
                flag
            ] = True

    return df_exp

# 上記2分類は背反である
# one-hot な表現を一つのカテゴリ変数にまとめる
# outside exon かつ上記の2分類に該当という場合もあり得る。その場合は outside exon を優先する
def assign_exon_category(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df["exon_category"] = "NAN"
    # 1. outside exon（最優先）
    df.loc[
        df["is_outside_common_exon_space"],
        "exon_category"
    ] = "outside_exon"
    # 2. outside でない exon のみ skipped 分類
    mask_inside = ~df["is_outside_common_exon_space"]
    df.loc[
        mask_inside & df["skipped_in_coding"],
        "exon_category"
    ] = "skipped_in_coding"
    df.loc[
        mask_inside & df["skipped_only_in_noncoding"],
        "exon_category"
    ] = "skipped_only_in_noncoding"
    return df

data_filtered = annotate_variant_counts(data_filtered)
df_exp = data_filtered.explode(["exontype", "exons", "frame", "exonlengths", "exon_position", "cds_info", "is_outside_common_exon_space"])

exon_df_alt_cds = extract_alternative_coding_exons(df_exp)
exon_df_alt_cds = annotate_exon_skipping(exon_df_alt_cds)
exon_df_alt_cds = assign_exon_category(exon_df_alt_cds)
exon_df_alt_cds = exon_df_alt_cds[exon_df_alt_cds["cds_info"] == "cds_exon"].drop_duplicates(subset=["chrom","exons"])
# %%
# exon_category ごとに in-frame 割合を計算する( outside exon, skipped_in_coding, skipped_only_in_noncoding )
inframe_summary = (
    exon_df_alt_cds
    .query("exon_category != 'NAN'")
    .groupby("exon_category")
    .agg(
        total_exons=("frame", "size"),
        inframe_exons=("frame", lambda x: (x == "in-frame").sum())
    )
    .assign(
        inframe_rate=lambda df: df["inframe_exons"] / df["total_exons"]
    )
    .reset_index()
)
# overall in-frame rate のを inframe_summary に追加
num_alternative_coding_exons = df_exp[
    (df_exp["exontype"].isin(["alternative", "unique-alternative"])) &
    (df_exp["cds_info"] == "cds_exon")
]["exons"].nunique()
num_inframe_alternative_coding_exons = exon_df_alt_cds[
    exon_df_alt_cds["frame"] == "in-frame"
]["exons"].nunique()
overall_inframe_rate = num_inframe_alternative_coding_exons / num_alternative_coding_exons
overall_summary = pd.DataFrame({
    "exon_category": ["overall_alternative_coding_exons"],
    "total_exons": [num_alternative_coding_exons],
    "inframe_exons": [num_inframe_alternative_coding_exons],
    "inframe_rate": [overall_inframe_rate]
})
inframe_summary = pd.concat([inframe_summary, overall_summary], ignore_index=True)
print(inframe_summary)

# %%
# inframe_summary を棒グラフ化
plot = (
    ggplot(inframe_summary, aes(x="exon_category", y="inframe_rate", fill="exon_category")) +
    geom_bar(stat="identity") +
    labs(title="% of in-frame exon by exon category in alternative coding exons", x="Exon Category", y="% of in-frame exon") +
    scale_fill_manual(values={
        "outside_exon": "#E69F00", #orange
        "skipped_in_coding": "#E69F00", #orange 
        "skipped_only_in_noncoding": "#E69F00", #orange
        "overall_alternative_coding_exons": "#009E73", #green
    }) +
    theme(
        axis_text_x=element_text(rotation=0, hjust=0.5, size=10),
        axis_title_y=element_text(size=10),
        axis_title_x=element_text(size=12),
        legend_title=element_text(size=10),
        legend_text=element_text(size=10),
        figure_size=(8,6)
    ) +
    scale_x_discrete(
        limits=[
                "skipped_only_in_noncoding",
                "skipped_in_coding",
                "outside_exon",
                "overall_alternative_coding_exons",
        ],
        labels={
            "overall_alternative_coding_exons": f"Overall alternative coding exons \n (n={num_alternative_coding_exons})",
            "outside_exon": f"Outside exon \n (n={inframe_summary.loc[inframe_summary['exon_category'] == 'outside_exon', 'total_exons'].values[0]})",
            "skipped_in_coding": f"Skipped in coding \n (n={inframe_summary.loc[inframe_summary['exon_category'] == 'skipped_in_coding', 'total_exons'].values[0]})",
            "skipped_only_in_noncoding": f"Constitutive in coding, \n skipped in non-coding \n (n={inframe_summary.loc[inframe_summary['exon_category'] == 'skipped_only_in_noncoding', 'total_exons'].values[0]})",
        }) +
    scale_y_continuous(limits=(0, 1), labels=lambda l: ["{:.0f}%".format(v * 100) for v in l]) +
    coord_flip()
)
plot.save("../data/mm39/inframe_ratio_by_exon_category_in_alternative_coding_exons.jpg", dpi=600)
display(plot)


# %%
skipped_in_coding_outframe = exon_df_alt_cds[(exon_df_alt_cds["exon_category"] == "skipped_in_coding") & (exon_df_alt_cds["frame"] == "out-frame")]
print(skipped_in_coding_outframe.tail(10))

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

labels = ['Genes have one exon', 'Genes have a single isoform', 'Genes have multiple isoforms']
sizes = [ genes_one_exon_count, genes_not_one_exon_but_one_variant_count, valid_genes]
colors = ['#D3D3D3', '#66CDAA', '#009E73']
plt.figure(figsize=(6, 6))
plt.pie(
    sizes, 
    labels=labels, 
    colors=colors, 
    autopct=lambda pct: autopct_with_count(pct, sizes), 
    startangle=90,
    textprops={'fontsize': '16'},
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
            "Alternative exon": "#CC79A7",  # ピンク
            "A3SS-long": "#CC79A7",  # オレンジ
            "A5SS-long": "#CC79A7",  # オレンジ
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