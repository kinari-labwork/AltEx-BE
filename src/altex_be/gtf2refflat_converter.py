from pathlib import Path
import logging
import re
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
def gtf_to_refflat(gtf_path: Path, output_path: Path, assembly_name: str) -> None:
    """
    GTFファイルをrefflat形式に変換して保存する関数
    GTFファイルのすべての行は、以下の構造になっている
    chrom  source  feature  start  end  score  strand  frame  attributes
    
    GTFの全体構造は以下のようになっている
    Gene A
        transcript 1 # トランスクリプトの場合は、transcript列が存在し、startとendがそれぞれトランスクリプトの開始位置と終了位置を示す
            exon 1
            exon 2
            CDS 1  #コーディングトランスクリプトの場合は、cds列が存在し、startとendがそれぞれCDSの開始位置と終了位置を示す
        transcript 2
            exon 1
            exon 2
            exon 3
    Gene B
        ...
    
    """
    output_refflat_path = output_path / f"converted_refflat_{assembly_name}.txt"
    if output_refflat_path.exists():
        logging.info(f"Converted refFlat file for {assembly_name} already exists. Skipping conversion.")
        return

    transcripts: dict[str, dict] = {}
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    gname_re = re.compile(r'gene_name "([^"]+)"')

    with open(gtf_path, "r", encoding="utf-8") as fh:
        for line in fh:
            # コメント行と空行をスキップ
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            # GTFの列数チェック
            if len(cols) < 9:
                continue
            # 列の展開
            chrom, source, feature, start_s, end_s, score, strand, frame, attr_str = cols
            start = int(start_s)
            end = int(end_s)
            tid_m = tid_re.search(attr_str)
            # transcript_idが見つからない場合はスキップ
            if not tid_m:
                continue
            tid = tid_m.group(1)
            # gene_nameの取得（存在しない場合は空文字）
            gname_m = gname_re.search(attr_str)
            gname = gname_m.group(1) if gname_m else ""

            # トランスクリプト情報の初期化または取得
            rec = transcripts.setdefault(
                tid,
                {
                    "gene_name": gname,
                    "chrom": chrom,
                    "strand": strand,
                    "tx_start": None,
                    "tx_end": None,
                    "cds_start": None,
                    "cds_end": None,
                    "exons": []
                }
            )

            # 優先して最初に見つかった gene_name/chrom/strand を保持
            if not rec["gene_name"] and gname:
                rec["gene_name"] = gname
            if not rec["chrom"]:
                rec["chrom"] = chrom
            if not rec["strand"]:
                rec["strand"] = strand

            if feature.lower() == "transcript":
                rec["tx_start"] = start if rec["tx_start"] is None else min(rec["tx_start"], start)
                rec["tx_end"] = end if rec["tx_end"] is None else max(rec["tx_end"], end)
            elif feature.lower() == "exon":
                rec["exons"].append((start - 1, end))  # UCSC: exonStarts 0-based
                rec["tx_start"] = start if rec["tx_start"] is None else min(rec["tx_start"], start)
                rec["tx_end"] = end if rec["tx_end"] is None else max(rec["tx_end"], end)
            elif feature.lower() == "cds":
                rec["cds_start"] = start if rec["cds_start"] is None else min(rec["cds_start"], start)
                rec["cds_end"] = end if rec["cds_end"] is None else max(rec["cds_end"], end)

    # 書き出し（トランスクリプトごとに1行）
    with open(output_refflat_path, "w", encoding="utf-8") as out:
        for tid, rec in transcripts.items():
            exons = sorted(rec["exons"], key=lambda x: x[0])
            exon_count = len(exons)
            exon_starts_s = ",".join(str(s) for s, e in exons) + ("," if exon_count else "")
            exon_ends_s = ",".join(str(e) for s, e in exons) + ("," if exon_count else "")

            chrom = rec["chrom"]
            if chrom and not str(chrom).startswith("chr"):
                chrom = f"chr{chrom}"

            tx_start_s = str(rec["tx_start"] - 1) if rec["tx_start"] is not None else ""
            tx_end_s = str(rec["tx_end"]) if rec["tx_end"] is not None else ""
            cds_start_s = str(rec["cds_start"] - 1) if rec["cds_start"] is not None else (tx_start_s if tx_start_s else "")
            cds_end_s = str(rec["cds_end"]) if rec["cds_end"] is not None else (tx_end_s if tx_end_s else "")

            out.write("\t".join([
                rec["gene_name"] or "",
                tid,
                chrom or "",
                rec["strand"] or "",
                tx_start_s,
                tx_end_s,
                cds_start_s,
                cds_end_s,
                str(exon_count),
                exon_starts_s,
                exon_ends_s
            ]) + "\n")
    return 
