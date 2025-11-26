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
data = pd.read_pickle("../data/target_exons_with_sgrna.pkl")

# %%
print(data.columns)

# %%
print(data.head())

# %%
matching_columns = list(data.filter(regex="acceptor_sgrna_target_sequence$").columns) + \
                   list(data.filter(regex="donor_sgrna_target_sequence$").columns)

designable_exons = data[
    data[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
    )
]

print(len(data))
print(len(designable_exons))
print(f'{len(designable_exons) / len(data) * 100} % of exons are designable')
print(designable_exons.head())

# %%
print("Matching columns:", matching_columns)

# %%
internal_exons_editable_both_sasd = data[(data["exon_position"] == "internal") & (data["exontype"].isin(["alternative", "unique-alternative"]))]
designable_internal_exons_editable_both_sasd = internal_exons_editable_both_sasd[
    internal_exons_editable_both_sasd[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
    )
]
print(len(internal_exons_editable_both_sasd))
print(len(designable_internal_exons_editable_both_sasd))
print(f'{len(designable_internal_exons_editable_both_sasd) / len(internal_exons_editable_both_sasd) * 100} % of exons are designable')

# %%
internal_exons_editable_only_sasd = data[(data["exon_position"] == "internal") & (data["exontype"].isin(["a3ss-long", "a5ss-long"]))]
designable_internal_exons_editable_only_sasd = internal_exons_editable_only_sasd[
    internal_exons_editable_only_sasd[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
    )
]
print(len(internal_exons_editable_only_sasd))
print(len(designable_internal_exons_editable_only_sasd))
print(f'{len(designable_internal_exons_editable_only_sasd) / len(internal_exons_editable_only_sasd) * 100} % of exons are designable')

# %%
first_exon_editable_sd = data[(data["exon_position"] == "first") & (data["exontype"].isin(["a5ss-long", "alternative", "unique-alternative"]))]
designable_first_exon_editable_sd = first_exon_editable_sd[
    first_exon_editable_sd[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
    )
]
print(len(first_exon_editable_sd))
print(len(designable_first_exon_editable_sd))
print(f'{len(designable_first_exon_editable_sd) / len(first_exon_editable_sd) * 100} % of exons are designable')

# %%
last_exon_editable_sa = data[(data["exon_position"] == "last") & (data["exontype"].isin(["a3ss-long", "alternative", "unique-alternative"]))]
designable_last_exon_editable_sa = last_exon_editable_sa[
    last_exon_editable_sa[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
    )
]
print(len(last_exon_editable_sa))
print(len(designable_last_exon_editable_sa))
print(f'{len(designable_last_exon_editable_sa) / len(last_exon_editable_sa) * 100} % of exons are designable')

# %%
first_exon_uneditable_sa = data[(data["exon_position"] == "first") & (data["exontype"].isin(["a3ss-long"]))]
designable_first_exon_uneditable_sa = first_exon_uneditable_sa[
    first_exon_uneditable_sa[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
    )
]
print(len(first_exon_uneditable_sa))
print(len(designable_first_exon_uneditable_sa))
print(f'{len(designable_first_exon_uneditable_sa) / len(first_exon_uneditable_sa) * 100} % of exons are designable')

# %%
last_exon_uneditable_sd = data[(data["exon_position"] == "last") & (data["exontype"].isin(["a5ss-long"]))]
designable_last_exon_uneditable_sd = last_exon_uneditable_sd[
    last_exon_uneditable_sd[matching_columns].apply(
        lambda row: any(isinstance(x, list) and len(x) > 0 for x in row), axis=1
    )
]
print(len(last_exon_uneditable_sd))
print(len(designable_last_exon_uneditable_sd))
print(f'{len(designable_last_exon_uneditable_sd) / len(last_exon_uneditable_sd) * 100} % of exons are designable')

# %%

# データフレームを作成
plot_data = pd.DataFrame({
    "Category": ["total target exon", "total target exon", "internal exon", "internal exon","internal a3ss or a5ss", "internal a3ss or a5ss", "first exon", "first exon", "last exon", "last exon"],
    "Designable": ["Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No"],
    "Count": [
        len(designable_exons),
        len(data) - len(designable_exons),
        len(designable_internal_exons_editable_both_sasd),
        len(internal_exons_editable_both_sasd) - len(designable_internal_exons_editable_both_sasd),
        len(internal_exons_editable_only_sasd),
        len(internal_exons_editable_only_sasd) - len(designable_internal_exons_editable_only_sasd),
        len(designable_first_exon_editable_sd),
        len(first_exon_editable_sd) - len(designable_first_exon_editable_sd),
        len(designable_last_exon_editable_sa),
        len(last_exon_editable_sa) - len(designable_last_exon_editable_sa),
    ]
})

plot_data["Category"] = pd.Categorical(
    plot_data["Category"],
    categories=[
        "total target exon", "internal exon", "internal a3ss or a5ss", "first exon", "last exon"
    ],
    ordered=True
)

total_counts = plot_data.groupby("Category")["Count"].transform("sum")
plot_data["Percentage"] = (plot_data["Count"] / total_counts * 100).round(1)

plot_data["y_position"] = plot_data.groupby("Category")["Count"].cumsum() - (plot_data["Count"] / 2)


# プロット
plot = (
    ggplot(plot_data, aes(x="Category", y="Count", fill="Designable")) +
    geom_bar(stat="identity", position="stack") +
    geom_text(
        aes(label="Percentage", y="Count"),  # y="Count" でバーの上に配置
        data=plot_data[plot_data["Designable"] == "Yes"],
        va="bottom",  # ラベルをバーの上に配置
        format_string="{:.1f}%", 
        size=12,
        color="black",
    ) +
    scale_fill_manual(
        values={
            "Yes": "#56B4E9",  # 青 (アクセシビリティ対応色)
            "No": "#E69F00"   # オレンジ (アクセシビリティ対応色)
        }
    ) +
    labs(
        title="Designable Exons by Category (target AID, CBE4max, ABE8e)",
        x="Exon Category",
        y="Number of Exons",
        fill="Designable"
    ) +
    theme(
        axis_text_x=element_text(size=11,rotation=45, hjust=1),
        figure_size=(8, 6),
    )
)

display(plot)  # Jupyter Notebookでの表示用

# %%
print(len(data))
print(len(designable_exons))
print(f'{len(designable_exons) / len(data) * 100} % of exons are designable')

# %%
def get_designable_exons_information(base_editor):
    # Get the relevant columns for the base editor
    acceptor_col = f"{base_editor}_acceptor_sgrna_target_sequence"
    donor_col = f"{base_editor}_donor_sgrna_target_sequence"
    # Mask: at least one of the two columns is a non-empty list
    mask = data[[acceptor_col, donor_col]].apply(
        lambda row: (isinstance(row[acceptor_col], list) and len(row[acceptor_col]) > 0) or
                    (isinstance(row[donor_col], list) and len(row[donor_col]) > 0),
        axis=1
    )
    base_editor_data = data[mask]
    result_statistics_df = pd.DataFrame({
        "Base_Editor": [base_editor],
        "Total_Exons": [len(data)],
        "Designable_Exons": [len(base_editor_data)],
        "Undesignable_Exons": [len(data) - len(base_editor_data)],
        "Percentage_of_Designable_Exons": [len(base_editor_data) / len(data) * 100]
    })
    return result_statistics_df

target_aid_statistics = get_designable_exons_information("target_aid")
target_be4max_statistics = get_designable_exons_information("be4max")
target_abe8e_statistics = get_designable_exons_information("abe8e")

total_exons = len(data)
designable_exons_count = len(designable_exons)
designable_percentage = (designable_exons_count / total_exons) * 100
overall_statistics = pd.DataFrame({
    "Base_Editor": ["Overall"],
    "Total_Exons": [total_exons],
    "Designable_Exons": [designable_exons_count],
    "Undesignable_Exons": [total_exons - designable_exons_count],
    "Percentage_of_Designable_Exons": [designable_percentage]
})

result_sammary_df = pd.concat([
    target_aid_statistics,
    target_be4max_statistics,
    target_abe8e_statistics,
    overall_statistics
], ignore_index=True)

print(result_sammary_df)



# %%
# データフレームを整形
plot_data = pd.DataFrame({
    "Base_Editor": ["target_aid", "target_aid", "be4max", "be4max", "abe8e", "abe8e", "Overall", "Overall"],
    "Designable": ["Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No"],
    "Count": [
        result_sammary_df.loc[0, "Designable_Exons"],
        result_sammary_df.loc[0, "Undesignable_Exons"],
        result_sammary_df.loc[1, "Designable_Exons"],
        result_sammary_df.loc[1, "Undesignable_Exons"],
        result_sammary_df.loc[2, "Designable_Exons"],
        result_sammary_df.loc[2, "Undesignable_Exons"],
        result_sammary_df.loc[3, "Designable_Exons"],
        result_sammary_df.loc[3, "Undesignable_Exons"],
    ]
})

# パーセンテージを計算
total_counts = plot_data.groupby("Base_Editor")["Count"].transform("sum")
plot_data["Percentage"] = (plot_data["Count"] / total_counts * 100).round(1)

# プロット
plot = (
    ggplot(plot_data, aes(x="Base_Editor", y="Count", fill="Designable")) +
    geom_bar(stat="identity", position="stack") +
    geom_text(
        aes(label="Percentage"),
        data=plot_data[plot_data["Designable"] == "Yes"],  # Designable のみラベルを表示
        position="stack",
        va="bottom",  # バーの上に配置
        size=10,
        color="black",
        format_string="{:.1f}%"  # パーセンテージを表示
    ) +
    scale_fill_manual(
        values={
            "Yes": "#56B4E9",  # 青
            "No": "#E69F00"   # オレンジ
        }
    ) +
    labs(
        title="Percentage of Designable Exons by Base Editor",
        x="Base Editor",
        y="Number of Target Exons",
        fill="Designable"
    ) +
    theme(
        axis_text_x=element_text(size=11, rotation=45, hjust=1),
        figure_size=(8, 6),
    )
)

display(plot)  # Jupyter Notebookでの表示用


