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



def modify_refFlat(refFlat: pd.DataFrame) -> pd.DataFrame: 
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
    # refflatのstartは0-baseでendは1-baseなので、毎回1を足す必要がない
    refFlat["exonlengths"] = refFlat.apply(
        lambda row: [end - start for start, end in zip(row["exonStarts"], row["exonEnds"])],
        axis=1)

    refFlat["exons"] = refFlat.apply(
        lambda row: list(zip(row["exonStarts"], row["exonEnds"])),
        axis=1
    )

    return refFlat

def classify_exon_type(
        target_exon: tuple[int:int], 
        all_transcripts: list[list[tuple[int,int]]], 
        ) -> str :
    """
    Purpose:
        タプル (start, end)の形式で与えられたexonが、ある遺伝子の全てのトランスクリプトに含まれるexonの(start, end)のタプルのリストに対して、
        どのようなsplicing eventに該当するかを判定する。
    Parameters:
        target_exon: タプル (start, end)
        all_transcripts: ある遺伝子の全てのトランスクリプトの (start, end)のタプルのリスト (次の関数で遺伝子ごとにグループ化してこの関数にinputする)
    Returns:
        exon_type: str
    """
    start, end = target_exon

    exact_match = 0
    start_match_only = False
    end_match_only = False
    overlap_without_startend_match = False
    

    for tx in all_transcripts:
        if target_exon in tx:
           exact_match += 1
           continue
        for ex in tx:
            if ex[0] == start and ex[1] != end: #ex[0]は比較対象のexonのstart,ex[1]は比較対象のexonのend
                start_match_only = True
            elif ex[1] == end and ex[0] != start: 
                end_match_only = True
            elif (ex[0] < end and ex[1] > start) and (ex[0] != start or ex[1] != end):
                overlap_without_startend_match = True

    total = len(all_transcripts)

    if exact_match == total: 
        return "constitutive" #そもそもsplicing variantがない場合は全てconstitutiveとなる
    elif exact_match > 1 and exact_match != total and not start_match_only and not end_match_only and not overlap_without_startend_match:
        return "cassette"
    elif start_match_only and not end_match_only:
        return "a3ss"
    elif end_match_only and not start_match_only:
        return "a5ss"
    elif overlap_without_startend_match:
        return "overlap"
    elif exact_match == 1 and exact_match != total:
        return "unique"
    else:
        return "other"
    
def classify_exons_per_gene(refflat: pd.DataFrame) -> pd.DataFrame:
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
            classify_exon_type(exon, all_transcripts)
            for exon in row["exons"]
        ]

    for gene, group in refflat.groupby("geneName"):
        group = group.copy()
        group["exontype"] = group.apply(
            lambda row: classify_row_exons(row, group), axis=1
        )
        result.append(group)

    return pd.concat(result, ignore_index=True)

