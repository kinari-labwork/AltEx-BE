import subprocess
from pathlib import Path
import pandas as pd
import re
import logging
from . import logging_config  # noqa: F401

def parse_attr(attr_str):
    """GTF attribute列を辞書にパース"""
    attrs = {}
    for kv in attr_str.strip().split(";"):
        if kv.strip():
            key, val = kv.strip().split(" ", 1)
            attrs[key] = val.strip('"')
    return attrs

# もはやこれだけでいい
def gtf_to_refflat(gtf_path, output_path):
    with open(gtf_path) as gtf, open(output_path, "w") as out:
        current_tx = None
        current_gene = None
        current_chrom = None
        current_strand = None
        exon_starts, exon_ends = [], []
        tx_start, tx_end = None, None
        cds_start, cds_end = None, None

        for line in gtf:
            if line.startswith("#") or not line.strip():
                continue
            # 一応いらない列も取得（使わないけど）
            chrom, source, feature, start, end, score, strand, frame, attr_str = line.strip().split("\t")
            start, end = int(start), int(end)
            attrs = parse_attr(attr_str)

            gene_symbol = attrs.get("gene_name")
            tx_id = attrs.get("transcript_id")

            # transcript が変わったら出力
            if tx_id != current_tx and current_tx is not None:
                exon_starts, exon_ends = zip(*sorted(zip(exon_starts, exon_ends)))
                out.write("\t".join([
                    current_gene,
                    current_tx,
                    current_chrom,
                    current_strand,
                    str(tx_start - 1),  # UCSC: 0-based start
                    str(tx_end),
                    str(cds_start - 1 if cds_start else tx_end), # CDSがない場合はCDS_start, endを両方tx_endにする
                    str(cds_end if cds_end else tx_end),
                    str(len(exon_starts)),
                    ",".join(map(str, exon_starts)) + ",",
                    ",".join(map(str, exon_ends)) + ","
                ]) + "\n")
                
                # パラメータを初期化する
                exon_starts, exon_ends = [], []
                tx_start, tx_end, cds_start, cds_end = None, None, None, None

            current_tx = tx_id
            current_gene = gene_symbol
            current_chrom = "chr" + str(chrom) # 染色体列に"chr"接頭辞を追加
            current_strand = strand

            if feature == "transcript":
                tx_start = start
                tx_end = end
            elif feature == "exon":
                exon_starts.append(start - 1)  # UCSC uses 0-based
                exon_ends.append(end)
            elif feature == "CDS":
                cds_start = start if cds_start is None else min(cds_start, start)
                cds_end = end if cds_end is None else max(cds_end, end)

        # 最後の transcript を出力
        if current_tx:
            exon_starts, exon_ends = zip(*sorted(zip(exon_starts, exon_ends)))
            out.write("\t".join([
                current_gene,
                current_tx,
                current_chrom,
                current_strand,
                str(tx_start - 1),
                str(tx_end),
                str(cds_start - 1 if cds_start else tx_start - 1),
                str(cds_end if cds_end else tx_end),
                str(len(exon_starts)),
                ",".join(map(str, exon_starts)) + ",",
                ",".join(map(str, exon_ends)) + ","
            ]) + "\n")
