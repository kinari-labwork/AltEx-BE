import pandas as pd
import logging
import argparse
from .. import logging_config # noqa: F401
from . import refflat_preprocessor

def process_refflat(refflat_path: str, interest_gene_list: list[str], parser: argparse.ArgumentParser) -> pd.DataFrame:
    """
    前処理から、興味のある遺伝子の抽出までを行う。
    """
    logging.info("-" * 50)
    logging.info("loading refFlat file...")
    refflat = pd.read_csv(
            refflat_path,
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
                "exonEnds",
            ],
        )
    
    logging.info("running processing of refFlat file...")
    refflat = refflat.drop_duplicates(subset=["name"], keep=False)
    refflat = refflat_preprocessor.preprocess_refflat(refflat, interest_gene_list)
    if refflat.empty :
        parser.error("No interest genes found in refFlat after preprocessing. Exiting...")
    # すべて constitutive exonでも設計対象とするが、exonが1つしかない遺伝子は対象外とする
    if not refflat_preprocessor.check_multiple_exon_existance(refflat, interest_gene_list) :
        parser.error("all of your interest genes are single-exon genes. AltEx-BE cannot process these genes. Exiting...")
    return refflat
