import pandas as pd

refflat = pd.read_csv(
    "data/refFlat.txt",
    sep="\t",
    header=None,
    names=[
        "geneName",
        "name",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds"])



def modify_refFlat(refFlat):
    """
    Purpose:
        exonStart exonEndが別々のカラムに格納されているので、(start, end)のタプルのリストに変換する。
    Parameters:
        refFlat: pd.DataFrame, refFlatのデータフレーム
    """
    # Convert the exonStarts and exonEnds columns to lists of integers
    refFlat["exonStarts"] = refFlat["exonStarts"].apply(
        lambda x: [int(i) for i in x.split(",") if i.strip() !='']) 
    refFlat["exonEnds"] = refFlat["exonEnds"].apply(
        lambda x: [int(i) for i in x.split(",") if i.strip() !=''])

    # Calculate the lengths of each exon
    refFlat["exonlengths"] = refFlat.apply(
        lambda row: [end - start for start, end in zip(row["exonStarts"], row["exonEnds"])],
        axis=1)

    refFlat["exons"] = refFlat.apply(
        lambda row: list(zip(row["exonStarts"], row["exonEnds"])),
        axis=1
    )

    return refFlat

def classify_exon_type(target_exon, all_transcripts, fuzzy_tolerance=5):
    """
    Purpose:
        タプル (start, end)の形式で与えられたexonが、ある遺伝子の全てのトランスクリプトに含まれるexonの(start, end)のタプルのリストに対して、
        どのようなsplicing eventに該当するかを判定する。
    Parameters:
        target_exon: タプル (start, end)
        all_transcripts: ある遺伝子の全てのトランスクリプトの (start, end)のタプルのリスト (次の関数で遺伝子ごとにグループ化してこの関数にinputする)
        fuzzy_tolerance: int, これ以下の塩基数のstart/endのずれをfuzzyとして規定する (デフォルトは5)
    Returns:
        exon_type: str
    """
    start, end = target_exon

    exact_match = 0
    start_match_only = 0
    end_match_only = 0
    fuzzy_match = 0
    flanked_exon_count = 0

    for tx in all_transcripts:
        if target_exon in tx:
            exact_match += 1
        else:
            for ex in tx:
                if ex[0] == start and ex[1] != end: #ex[0]は比較対象のexonのstart,ex[1]は比較対象のexonのend
                    start_match_only += 1
                elif ex[1] == end and ex[0] != start: 
                    end_match_only += 1
                elif abs(ex[0] - start) <= fuzzy_tolerance and abs(ex[1] - end) <= fuzzy_tolerance: #fuzzy_toleranceパラメータで設定した容認する基準以下のずれが5',3'両方にあるexonを判定する
                    fuzzy_match += 1

    total = len(all_transcripts)

    if exact_match == total: 
        return "constitutive" #そもそもsplicing variantがない場合は全てconstitutiveとなる
    elif exact_match > 1 and exact_match != total:
        return "cassette"
    elif start_match_only > 0 and end_match_only == 0:
        return "a5ss"
    elif end_match_only > 0 and start_match_only == 0:
        return "a3ss"
    elif fuzzy_match > 0:
        return "fuzzy"
    elif exact_match == 1 and exact_match != total:
        return "unique"
    else:
        return "other"
    
def classify_exons_per_gene(refflat, fuzzy_tolerance=3):
    refflat = refflat.copy()
    result = []

    """
    classify_exon_types()は入力が(start, end)のタプルであることを前提としているので、
    refflatのデータフレームをそのまま入力しても正しく動作しない。
    refflatの各行に対して、exons列を(start, end)のタプルのリストに変換する。
    そしてclassify_row_exons()関数を適用して、各行のexons列に対して分類を行う。
    その結果を新しい列"exontype"に追加する。
    そのリストをpd.concat()で結合して、最終的なDataFrameを返す。
    """
    def classify_row_exons(row, gene_group):
        all_transcripts = gene_group["exons"].tolist()
        return [
            classify_exon_type(exon, all_transcripts, fuzzy_tolerance)
            for exon in row["exons"]
        ]

    for gene, group in refflat.groupby("geneName"):
        group = group.copy()
        group["exontype"] = group.apply(
            lambda row: classify_row_exons(row, group), axis=1
        )
        result.append(group)

    return pd.concat(result, ignore_index=True)

