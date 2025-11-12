import subprocess
from pathlib import Path
import pandas as pd
import re
import logging
from . import logging_config  # noqa: F401

def run_gtf2genepred(gtf_path: Path, output_path: Path):
    """
    Convert GTF file to genePred format using gtfToGenePred tool.
    Args:
        gtf_path (Path): Path to the input GTF file.
        output_path (Path): Path to the output genePred file.
    returns: None
    """
    cmd = ["gtfToGenePred", str(gtf_path), str(output_path)]
    logging.info(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    logging.info(f"Converted GTF to genePred: {output_path}")

def convert_gtf_to_refflat(gtf_path: Path, output_path: Path):
    """
    Convert GTF file to refFlat format.
    Args:
        gtf_path (Path): Path to the input GTF file.
        refflat_path (Path): Path to the output refFlat file.
    returns: pd.DataFrame the refFlat dataframe but without geneName column
    """
    genepred_path = output_path.with_suffix('.genePred')
    genepred = pd.read_csv(genepred_path, sep="\t", header=None)
    refflat_without_genesymbol = genepred.iloc[:, [12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]]
    # 染色体列の"chr"接頭辞を追加
    refflat_without_genesymbol[1] = "chr" + refflat_without_genesymbol[1].astype(str)
    return refflat_without_genesymbol

def add_genesymbol_to_imcomplete_refflat(refflat_without_genesymbol: pd.DataFrame, gtf_path: Path) -> pd.DataFrame:
    """
    GTFファイルから遺伝子記号を取得し、refFlatデータフレームに追加する。

    Args:
        refflat (pd.DataFrame): refFlatデータフレーム（geneName列が欠損している場合がある）
        gtf_path (Path): GTFファイルのパス
    Returns:
        pd.DataFrame: geneName列が右端に追加されたrefFlatデータフレーム
    """
    dict_transcript_id_to_gene = {}
    with gtf_path.open("r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "transcript":
                continue
            attribute_field = fields[8]
            transcript_id_match = re.search(r'transcript_id "([^"]+)"', attribute_field)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attribute_field)
            if transcript_id_match and gene_name_match:
                transcript_id = transcript_id_match.group(1)
                gene_name = gene_name_match.group(1)
                dict_transcript_id_to_gene[transcript_id] = gene_name

    map_df = pd.DataFrame({
        "transcript_id": list(dict_transcript_id_to_gene.keys()),
        "geneName": list(dict_transcript_id_to_gene.values())
    })

    transcript_col = refflat_without_genesymbol.columns[0]
    merged = refflat_without_genesymbol.merge(map_df, left_on=transcript_col, right_on="transcript_id", how="left")
    merged = merged.iloc[:, 13, :]  # geneName列を右端に移動
    return merged
