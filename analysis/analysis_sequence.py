import pandas as pd
from plotnine import *

pd.set_option('display.width', 200)

data = pd.read_pickle("./data/target_exons_with_acceptor_and_donor_sequence.pkl")

print(data.head())

# エキソンタイプの仕分け
internal_plus_exons = (data["exon_position"] == "internal") & (data["strand"] == "+")
internal_minus_exons = (data["exon_position"] == "internal") & (data["strand"] == "-")
internal_exons = internal_plus_exons | internal_minus_exons
first_plus_exons = (data["exon_position"] == "first") & (data["strand"] == "+")
first_minus_exons = (data["exon_position"] == "first") & (data["strand"] == "-")
first_exons = first_plus_exons | first_minus_exons
last_plus_exons = (data["exon_position"] == "last") & (data["strand"] == "+")
last_minus_exons = (data["exon_position"] == "last") & (data["strand"] == "-")
last_exons = last_plus_exons | last_minus_exons

# exonpositionの数の確認
print("\n--- general information ---\n")
print("Number of exons in the data:", len(data))
print("Number of internal exons:", internal_exons.sum())
print("Number of first exons:", first_exons.sum())
print("Number of last exons:", last_exons.sum())

# 理論的にsplicingが起こるはずのjunctionについてGT-AGルールの確認 (firstはacceptorが存在する必要がなく、lastはdonorが存在する必要がない)
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

#　ここからはforループで処理する
exon_types = [
    ("internal plus", internal_plus_exons),
    ("internal minus", internal_minus_exons),
    ("first plus", first_plus_exons),
    ("first minus", first_minus_exons),
    ("last plus", last_plus_exons),
    ("last minus", last_minus_exons),
]

for label, mask in exon_types:
    acc_seq = data.loc[mask, "acceptor_sequence"]
    don_seq = data.loc[mask, "donor_sequence"]
    is_AG = acc_seq.str[23:25] == "AG"
    is_ag = acc_seq.str[23:25] == "ag"
    is_GT = don_seq.str[25:27] == "GT"
    is_gt = don_seq.str[25:27] == "gt"
    total = mask.sum()
    print(f"\n--- {label} exons ---\n")
    print(f"Number of {label} exons: {total}")
    print(f"Number with AG acceptor sequence: {is_AG.sum()}")
    print(f"Number with GT donor sequence: {is_GT.sum()}")
    print(f"Number with ag acceptor sequence: {is_ag.sum()}")
    print(f"Number with gt donor sequence: {is_gt.sum()}")
    print(f"Number with AG or ag acceptor sequence: {(is_AG | is_ag).sum()}")
    print(f"Number with GT or gt donor sequence: {(is_GT | is_gt).sum()}")
    print(f"Percentage with AG acceptor sequence: {((is_AG | is_ag).sum() / total * 100)}")
    print(f"Percentage with GT donor sequence: {((is_GT | is_gt).sum() / total * 100)}")