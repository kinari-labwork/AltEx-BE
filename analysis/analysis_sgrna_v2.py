# %%
import pandas as pd
import numpy as np
from plotnine import *

# %%
pd.set_option('display.width', 100000000)
pd.set_option("display.max_columns", None)  # すべての列を表示
pd.set_option("display.max_colwidth", None)
P = print

# %%
def create_df_for_plotting(assembly_input_dir, base_editors):
    sgrna_df = pd.read_pickle(f"{assembly_input_dir}/sgrna_designed_all_genes.pkl")
    matching_columns = list(sgrna_df.filter(regex=f"{base_editors}_acceptor_sgrna_target_sequence$").columns) + list(sgrna_df.filter(regex=f"{base_editors}_donor_sgrna_target_sequence$").columns)

    total_exons = len(sgrna_df)

    constitutive_exons = sgrna_df[(sgrna_df["exontype"].isin(["constitutive"]))]
    designable_constitutive_exons = len(
        constitutive_exons[
        constitutive_exons[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
        )
    ])

    alternative_or_unique_alternative_exons = sgrna_df[(sgrna_df["exontype"].isin(["alternative", "unique-alternative"]))]
    designable_alternative_or_unique_alternative_exons = len(
        alternative_or_unique_alternative_exons[
        alternative_or_unique_alternative_exons[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
        )
    ])

    a5ss_long_exons = sgrna_df[(sgrna_df["exontype"].isin(["a5ss-long"]))]
    designable_a5ss_long_exons = len(
        a5ss_long_exons[
        a5ss_long_exons[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
        )
    ])

    a3ss_long_exons = sgrna_df[(sgrna_df["exontype"].isin(["a3ss-long"]))]
    designable_a3ss_long_exons = len(
        a3ss_long_exons[
        a3ss_long_exons[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
        )
    ])


    plot_data = pd.DataFrame({
    "category": [f"{base_editors}_constitutive_exons", f"{base_editors}_alternative_or_unique_alternative_exons", f"{base_editors}_a5ss_long_exons", f"{base_editors}_a3ss_long_exons"], 
    "base_editor": [base_editors] * 4,
    "total_targetable_exons": [total_exons] * 4,
    "total_exons_per_splicing_events": [
        len(constitutive_exons),
        len(alternative_or_unique_alternative_exons),
        len(a5ss_long_exons),
        len(a3ss_long_exons),
    ],
    "designable_count": [
        designable_constitutive_exons,
        designable_alternative_or_unique_alternative_exons,
        designable_a5ss_long_exons,
        designable_a3ss_long_exons,
    ],
    "designable_percentage": [
        designable_constitutive_exons / len(constitutive_exons) * 100 if len(constitutive_exons) > 0 else 0,
        designable_alternative_or_unique_alternative_exons / len(alternative_or_unique_alternative_exons) * 100 if len(alternative_or_unique_alternative_exons) > 0 else 0,
        designable_a5ss_long_exons / len(a5ss_long_exons) * 100 if len(a5ss_long_exons) > 0 else 0,
        designable_a3ss_long_exons / len(a3ss_long_exons) * 100 if len(a3ss_long_exons) > 0 else 0,
    ],
    })
    plot_data["category"] = pd.Categorical(
    plot_data["category"],
    categories=[
        f"{base_editors}_constitutive_exons", f"{base_editors}_alternative_or_unique_alternative_exons", f"{base_editors}_a5ss_long_exons", f"{base_editors}_a3ss_long_exons"
    ],
    ordered=True
    )

    return plot_data

# %%
assembly_list = ["hg38", "mm39", "rn7", "danRer11"]
base_editors = ["target_aid_ngg", "be4max_ngg", "abe8e_ngg", "target_aid_ng", "be4max_ng", "abe8e_ng"]

# %%
plot_data = pd.DataFrame()
for assembly_name in assembly_list:   
    for base_editor in base_editors:
        assembly_input_dir = f"../data/{assembly_name}"
        plot_data_sub = create_df_for_plotting(assembly_input_dir, base_editor)
        plot_data_sub["assembly"] = assembly_name
        plot_data = pd.concat([plot_data, plot_data_sub], ignore_index=True)

# %%
print(plot_data)

# %%
plot_data_mm39 = plot_data[plot_data["assembly"] == "mm39"]

# %%
empty_rows = pd.DataFrame({
    "category": [" ", "  ", "   ", "    ", "     "],
    "base_editor": [None] * 5,
    "total_targetable_exons": [None] * 5,
    "total_exons_per_splicing_events": [None] * 5,
    "designable_count": [None] * 5,
    "designable_percentage": [None] * 5,
})

plot_data_mm39 = pd.concat([plot_data_mm39, empty_rows], ignore_index=True)

plot = (
    ggplot(plot_data_mm39, aes(x="category", y="designable_percentage", fill="base_editor")) +
    geom_bar(stat="identity", position="stack") +
    labs(
        title=f"Designable exons in mm39 across Base Editors and Splicing events",
        x="Exon Category",
        y="Percentage of Designable Exons (%)",
        fill="Base Editor"
    ) +
    theme(
        axis_text_x=element_text(size=10, rotation=90, hjust=0.5),
        figure_size=(8, 6),
    ) +
    scale_fill_manual(
        values={
            "target_aid_ngg": "#56B4E9",  # 青
            "target_aid_ng": "#E69F00",  # オレンジ
            "be4max_ngg": "#009E73",   # 緑
            "be4max_ng": "#CC79A7",  # ピンク
            "abe8e_ngg": "#D55E00",  # 赤
            "abe8e_ng": "#575750FF",  # グレー
        }
    ) +
    
    scale_x_discrete(
        limits =["target_aid_ngg_constitutive_exons", "target_aid_ngg_alternative_or_unique_alternative_exons", "target_aid_ngg_a5ss_long_exons", "target_aid_ngg_a3ss_long_exons",
                " ",
                "target_aid_ng_constitutive_exons", "target_aid_ng_alternative_or_unique_alternative_exons", "target_aid_ng_a5ss_long_exons", "target_aid_ng_a3ss_long_exons",
                "  ",
                "be4max_ngg_constitutive_exons", "be4max_ngg_alternative_or_unique_alternative_exons", "be4max_ngg_a5ss_long_exons", "be4max_ngg_a3ss_long_exons",
                "   ",
                "be4max_ng_constitutive_exons", "be4max_ng_alternative_or_unique_alternative_exons", "be4max_ng_a5ss_long_exons", "be4max_ng_a3ss_long_exons",
                "    ",
                "abe8e_ngg_constitutive_exons", "abe8e_ngg_alternative_or_unique_alternative_exons", "abe8e_ngg_a5ss_long_exons", "abe8e_ngg_a3ss_long_exons",
                "     ",
                "abe8e_ng_constitutive_exons", "abe8e_ng_alternative_or_unique_alternative_exons", "abe8e_ng_a5ss_long_exons", "abe8e_ng_a3ss_long_exons"
        ],
        labels={
            "target_aid_ngg_constitutive_exons": "Constitutive Exons",
            "target_aid_ngg_alternative_or_unique_alternative_exons": "Alternative Exons",
            "target_aid_ngg_a5ss_long_exons": "A5SS-Long Exons",
            "target_aid_ngg_a3ss_long_exons": "A3SS-Long Exons",
            "target_aid_ng_constitutive_exons": "Constitutive Exons",
            "target_aid_ng_alternative_or_unique_alternative_exons": "Alternative Exons",
            "target_aid_ng_a5ss_long_exons": "A5SS-Long Exons",
            "target_aid_ng_a3ss_long_exons": "A3SS-Long Exons",
            "be4max_ngg_constitutive_exons": "Constitutive Exons",
            "be4max_ngg_alternative_or_unique_alternative_exons": "Alternative Exons",
            "be4max_ngg_a5ss_long_exons": "A5SS-Long Exons",
            "be4max_ngg_a3ss_long_exons": "A3SS-Long Exons",
            "be4max_ng_constitutive_exons": "Constitutive Exons",
            "be4max_ng_alternative_or_unique_alternative_exons": "Alternative Exons",
            "be4max_ng_a5ss_long_exons": "A5SS-Long Exons",
            "be4max_ng_a3ss_long_exons": "A3SS-Long Exons",
            "abe8e_ngg_constitutive_exons": "Constitutive Exons",
            "abe8e_ngg_alternative_or_unique_alternative_exons": "Alternative Exons",
            "abe8e_ngg_a5ss_long_exons": "A5SS-Long Exons",
            "abe8e_ngg_a3ss_long_exons": "A3SS-Long Exons",
            "abe8e_ng_constitutive_exons": "Constitutive Exons",
            "abe8e_ng_alternative_or_unique_alternative_exons": "Alternative Exons",
            "abe8e_ng_a5ss_long_exons": "A5SS-Long Exons",
            "abe8e_ng_a3ss_long_exons": "A3SS-Long Exons",
        }
    )
)
display(plot)

# %%
## アセンブリごとに、デザイン可能なすべてのエキソンを足す
all_gene_score = plot_data.groupby(["assembly", "base_editor", "total_targetable_exons"])["designable_count"].sum().reset_index()
all_gene_score["category"] = all_gene_score["assembly"] + "-" + all_gene_score["base_editor"].str.replace("_", " ")
all_gene_score["designable_percentage"] = all_gene_score["designable_count"] / all_gene_score["total_targetable_exons"] * 100

# %%
empty_rows = pd.DataFrame({
    "Category": [" ", "  ", "   "],  # 空白用のダミー値
    "Percentage_Designable": [None, None, None],  # 値は None
    "Assembly": [None, None, None]  # 他の列もダミー値
})
all_gene_score = pd.concat([all_gene_score, empty_rows], ignore_index=True)

plot = (
    ggplot(all_gene_score, aes(x="category", y="designable_percentage", fill="assembly")) +
    geom_bar(stat="identity", position="stack") +
    labs(
        title=f"Designable alternative-spliced exons across species and Base Editors",
        x="Exon Category",
        y="Percentage of Designable Exons (%)",
        fill="Assembly"
    ) +
    theme(
        axis_text_x=element_text(size=11, rotation=90, hjust=0.5),
        figure_size=(8, 6),
    ) +
    scale_fill_manual(
        values={
            "hg38": "#56B4E9",  # 青
            "mm39": "#E69F00",  # オレンジ
            "rn7": "#009E73",   # 緑
            "danRer11": "#CC79A7"  # ピンク
        }
    ) +
    scale_x_discrete(
        limits=[
            "hg38-target aid ngg", "hg38-target aid ng", "hg38-be4max ngg", "hg38-be4max ng", "hg38-abe8e ngg", "hg38-abe8e ng",
            "  ",
            "mm39-target aid ngg", "mm39-target aid ng", "mm39-be4max ngg", "mm39-be4max ng", "mm39-abe8e ngg", "mm39-abe8e ng",
            "   ",
            "rn7-target aid ngg", "rn7-target aid ng", "rn7-be4max ngg", "rn7-be4max ng", "rn7-abe8e ngg", "rn7-abe8e ng",
            "    ",
            "danRer11-target aid ngg", "danRer11-target aid ng", "danRer11-be4max ngg", "danRer11-be4max ng", "danRer11-abe8e ngg", "danRer11-abe8e ng"
        ],
        
        labels={
            "hg38-target aid ngg": "Target-AID NGG",
            "hg38-target aid ng": "Target-AID NG",
            "hg38-be4max ngg": "BE4max NGG",
            "hg38-be4max ng": "BE4max NG",
            "hg38-abe8e ngg": "ABE8e NGG",
            "hg38-abe8e ng": "ABE8e NG",
            "mm39-target aid ngg": "Target-AID NGG",
            "mm39-target aid ng": "Target-AID NG",
            "mm39-be4max ngg": "BE4max NGG",
            "mm39-be4max ng": "BE4max NG",
            "mm39-abe8e ngg": "ABE8e NGG",
            "mm39-abe8e ng": "ABE8e NG",
            "rn7-target aid ngg": "Target-AID NGG",
            "rn7-target aid ng": "Target-AID NG",
            "rn7-be4max ngg": "BE4max NGG",
            "rn7-be4max ng": "BE4max NG",
            "rn7-abe8e ngg": "ABE8e NGG",
            "rn7-abe8e ng": "ABE8e NG",
            "danRer11-target aid ngg": "Target-AID NGG",
            "danRer11-target aid ng": "Target-AID NG",
            "danRer11-be4max ngg": "BE4max NGG",
            "danRer11-be4max ng": "BE4max NG",
            "danRer11-abe8e ngg": "ABE8e NGG",
            "danRer11-abe8e ng": "ABE8e NG"
        }
    )
)

display(plot)
plot.save("../data/designable_exons_across_species_and_base_editors_colored.png", dpi=300)


