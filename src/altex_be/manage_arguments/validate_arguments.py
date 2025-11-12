from ..class_def.base_editors import BaseEditor
import argparse
import logging
from pathlib import Path
from .. import logging_config  # noqa: F401

def is_input_output_directories(
    refflat_path: Path, 
    fasta_path: Path, 
    output_directory: Path, 
    parser: argparse.ArgumentParser
) -> None:
    if not refflat_path.is_file():
        parser.error(f"The provided refFlat file '{refflat_path}' does not exist.")
    if not fasta_path.is_file():
        parser.error(f"The provided FASTA file '{fasta_path}' does not exist.")
    if not output_directory.is_dir():
        parser.error(f"The provided output directory '{output_directory}' does not exist.")
    return

def is_base_editors_provided(base_editors: dict[str, BaseEditor], parser: argparse.ArgumentParser) -> None:
    if not base_editors:
        parser.error("No base editors specified. Please provide at least one base editor.")

def is_interest_genes_provided(interest_gene_list: list[str], parser: argparse.ArgumentParser) -> None:
    if not interest_gene_list:
        parser.error("Please provide at least one interest gene symbol or Refseq ID.")

def load_supported_assemblies() -> list[str]:
    """
    パッケージ内のcrispr_direct_supported_assemblies.txtを読み込み、アセンブリ名リストを返す
    """
    # このファイルと同じディレクトリにあるtxtを参照
    txt_path = Path(__file__).parent / "crispr_direct_supported_assemblies.txt"
    with open(txt_path, encoding="utf-8") as f:
        supported_assemblies = {line.strip() for line in f if line.strip() and not line.startswith("#")}
    return supported_assemblies

def is_supported_assembly_name_in_crispr_direct(assembly_name:str) -> None:
    """
    Purpose: ユーザーが入力したアセンブリ名がCRISPRdirectでサポートされているかを確認する。
    Parameter: 
            assembly_name (str): ユーザーが入力したアセンブリ名
            supported_assemblies (list[str]): CRISPRdirectでサポートされているアセンブリ名のリスト 長いので外部txtとして保存
    return : bool
    """
    supported_assemblies = load_supported_assemblies()
    if assembly_name not in supported_assemblies: 
        return logging.warning(f"your_assembly : {assembly_name} is not supported by CRISPRdirect. please see <https://crispr.dbcls.jp/doc/>")
    return None

def validate_arguments(
    refflat_path: Path, 
    fasta_path: Path, 
    output_directory: Path, 
    interest_gene_list: list[str], 
    base_editors: dict[str, BaseEditor], 
    assembly_name: str, 
    parser: argparse.ArgumentParser
) -> None:
    """
    引数の妥当性を検証するラッパー関数
    ここで検証が通らなかった場合、parser.errorで終了する
    """
    is_input_output_directories(refflat_path, fasta_path, output_directory, parser)
    is_base_editors_provided(base_editors, parser)
    is_interest_genes_provided(interest_gene_list, parser)
    is_supported_assembly_name_in_crispr_direct(assembly_name)
    return None