from __future__ import (
    annotations,  # python 3.8以下の型ヒントの頭文字は大文字でないといけない
)

import pandas as pd


def classify_splicing_event(
    target_exon: tuple[int:int],
    all_transcripts: list[list[tuple[int, int]]],
) -> str:
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
    target_start, target_end = target_exon

    exact_match = 0
    has_start_match_only = False
    has_start_match_only_and_exist_later_end = False
    has_end_match_only = False
    has_end_match_only_and_exist_earlier_start = False
    has_overlap_without_startend_match = False

    for transcript in all_transcripts:
        if target_exon in transcript:
            exact_match += 1
            continue
        for exon in transcript:
            if (
                exon[0] == target_start and exon[1] != target_end
            ):  # ex[0]は比較対象のexonのstart,ex[1]は比較対象のexonのend
                has_start_match_only = (
                    True  # start だけ他のエキソンと一致し、 end は一致しない
                )
                if exon[1] > target_end: #startが一致するが、endが他のエキソンのendよりも小さい場合
                    has_start_match_only_and_exist_later_end = True
            elif exon[1] == target_end and exon[0] != target_start:
                has_end_match_only = (
                    True  # endだけ他のエキソンと一致し、startは一致しない場合
                )
                if exon[0] < target_start: # endが一致するが、startが他のエキソンのstartよりも大きい場合
                    has_end_match_only_and_exist_earlier_start = True
            elif (exon[0] < target_end and exon[1] > target_start) and (exon[0] != target_start or exon[1] != target_end):
                has_overlap_without_startend_match = True  # start, endどちらも他のエキソンと一致しないが、他のエキソンと1塩基以上の重複が生じている

    total = len(all_transcripts)

    if exact_match == total:
        return (
            "constitutive"  # そもそもsplicing variantがない場合は全てconstitutiveとなる
        )
    elif (
        exact_match > 1
        and exact_match != total
        and not has_start_match_only
        and not has_end_match_only
        and not has_overlap_without_startend_match
    ):
        return "alternative"  # 2つ以上のトランスクリプトに存在するが、全ての転写物には存在しないエキソン
    if has_start_match_only and not has_end_match_only and has_start_match_only_and_exist_later_end:
        return "a5ss-short" # ほかのエキソンとstartが一致するが、endはstartが一致する他のエキソンのendよりも小さい場合 例: [100,200],[100,300]の場合、[100,200]がa5ss-short, [100,300]がa5ss-long
    if has_start_match_only and not has_end_match_only and not has_start_match_only_and_exist_later_end:
        return "a5ss-long" # ほかのエキソンとstartが一致するが、endはstartが一致するエキソンの中で一番長い
    if has_end_match_only and not has_start_match_only and has_end_match_only_and_exist_earlier_start:
        return "a3ss-short" # ほかのエキソンとendが一致するが、startはendが一致する他のエキソンのstartよりも大きい場合 例: [100,200],[150,200]の場合、[100,200]がa3ss-long, [150,200]がa3ss-short
    if has_end_match_only and not has_start_match_only and not has_end_match_only_and_exist_earlier_start:
        return "a3ss-long"
    if has_start_match_only and has_end_match_only:
        return "intron_retention"
    if has_overlap_without_startend_match:
        return "overlap"
    if exact_match == 1 and exact_match != total:
        return "unique-alternative"  # 他のトランスクリプトには全く見られないエキソン
    else:
        return "other"
    

def classify_splicing_events_per_gene(refflat: pd.DataFrame) -> pd.DataFrame:
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
        return [classify_splicing_event(exon, all_transcripts) for exon in row["exons"]]

    for gene, group in refflat.groupby("geneName"):
        group = group.copy()
        group["exontype"] = group.apply(
            lambda row: classify_row_exons(row, group), axis=1
        )
        result.append(group)

    return pd.concat(result, ignore_index=True)


def flip_a3ss_a5ss_on_minus_strand(classified_refflat: pd.DataFrame) -> pd.DataFrame:
    """
    purpose:
        strandが-の遺伝子では、遺伝子の方向性が逆転するため、転写物のa3ssとa5ssが入れ替わるが、
        classify_exons_per_gene()では方向性を考慮していない。したがってstrandが-のa3ssとa5ssを入れ替える
    Parameter:
        classified_refflat: exontype列とstrand列を持つpd.DataFrame
    Return
        pd.DataFrame
    """

    flip_dict = {"a3ss-short": "a5ss-short", "a5ss-short": "a3ss-short", "a3ss-long": "a5ss-long", "a5ss-long": "a3ss-long"}
    mask = classified_refflat["strand"] == "-"

    # apply: リスト内の a3ss/a5ss を flip_dict で置換
    classified_refflat.loc[mask, "exontype"] = classified_refflat.loc[
        mask, "exontype"
    ].apply(lambda types: [flip_dict.get(t, t) for t in types])

    return classified_refflat

def classify_splicing_events(refflat: pd.DataFrame) -> pd.DataFrame:
    """
    このモジュールのwrap関数
    """
    classified_refflat = classify_splicing_events_per_gene(refflat)
    classified_refflat = flip_a3ss_a5ss_on_minus_strand(classified_refflat)
    return classified_refflat
