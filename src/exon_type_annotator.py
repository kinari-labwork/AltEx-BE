from __future__ import annotations #python 3.8以下の型ヒントの頭文字は大文字でないといけない
import pandas as pd


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
                start_match_only = True #start だけ他のエキソンと一致し、 end は一致しない
            elif ex[1] == end and ex[0] != start: 
                end_match_only = True #endだけ他のエキソンと一致し、startは一致しない場合
            elif (ex[0] < end and ex[1] > start) and (ex[0] != start or ex[1] != end):
                overlap_without_startend_match = True #start, endどちらも他のエキソンと一致しないが、他のエキソンと1塩基以上の重複が生じている

    total = len(all_transcripts)

    if exact_match == total: 
        return "constitutive" #そもそもsplicing variantがない場合は全てconstitutiveとなる
    elif exact_match > 1 and exact_match != total and not start_match_only and not end_match_only and not overlap_without_startend_match: 
        return "skipped" #2つ以上のトランスクリプトに存在するが、全ての転写物には存在しないエキソン
    elif start_match_only and not end_match_only:
        return "a5ss"
    elif end_match_only and not start_match_only:
        return "a3ss"
    elif start_match_only and end_match_only:
        return "intron_retention" 
    elif overlap_without_startend_match:
        return "overlap"
    elif exact_match == 1 and exact_match != total:
        return "unique" #他のトランスクリプトには全く見られないエキソン
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

def flip_a3ss_a5ss_in_minus_strand(classified_refflat: pd.DataFrame) -> pd.DataFrame:
    """
    purpose:
        strandが-の遺伝子では、遺伝子の方向性が逆転するため、転写物のa3ssとa5ssが入れ替わるが、
        classify_exons_per_gene()では方向性を考慮していない。したがってstrandが-のものだけ入れ替える
    Parameter:
        classified_refflat: exontype列とstrand列を持つpd.DataFrame
    Return
        pd.DataFrame
    """

    flip_dict = {"a3ss": "a5ss", "a5ss": "a3ss"}
    classified_refflat = classified_refflat.copy()
    mask = classified_refflat["strand"] == "-"
    
    # apply: リスト内の a3ss/a5ss を flip_dict で置換
    classified_refflat.loc[mask, "exontype"] = classified_refflat.loc[mask, "exontype"].apply(
        lambda types: [flip_dict.get(t, t) for t in types]
    )

    return classified_refflat