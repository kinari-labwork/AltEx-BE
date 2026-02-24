import argparse
import pandas as pd
import logging
from pathlib import Path
from .. class_def.base_editors import BaseEditor, PRESET_BASE_EDITORS
from .. import logging_config  # noqa: F401

def parse_gene_file(gene_file: Path) -> list[str] | None:
    if not gene_file:
        return None
    if gene_file.suffix.lower() not in [".txt", ".tsv", ".csv"]:
        raise ValueError("Unsupported file extension for gene file. Use .txt, .tsv, or .csv")
    with open(gene_file, "r") as f:
        genes_from_file = [line.strip() for line in f if line.strip()] #空の行は if line.strip がFalseになるので除外できる
    return genes_from_file

def parse_path_from_args(args: argparse.Namespace):
    refflat_path = Path(args.refflat_path)  if args.refflat_path else None
    gtf_path = Path(args.gtf_path) if args.gtf_path else None
    fasta_path = Path(args.fasta_path)
    output_directory = Path(args.output_dir)
    return refflat_path, gtf_path, fasta_path, output_directory

def parse_assembly_name_from_args(args: argparse.Namespace) -> str:
    return str(args.assembly_name)

def parse_genes_from_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> list[str]:
    genes_from_file = parse_gene_file(Path(args.gene_file)) if args.gene_file else []
    gene_symbols = args.gene_symbols if args.gene_symbols is not None else []
    refseq_ids = args.refseq_ids if args.refseq_ids is not None else []
    ensembl_ids = args.ensembl_ids if args.ensembl_ids is not None else []
    interest_gene_list = gene_symbols + refseq_ids + ensembl_ids + genes_from_file
    return interest_gene_list

def parse_base_editors_from_file(
    args: argparse.Namespace, 
    parser: argparse.ArgumentParser, 
    ) -> dict[str, BaseEditor] | None:
    """
    base editorの情報を含むファイルのパスを示す引数を受け取り、BaseEditorのリストを返す。
    csvまたはtxt, tsv形式のファイルをサポートする
    """
    expected_columns = [
    "base_editor_name",
    "pam_sequence",
    "editing_window_start",
    "editing_window_end",
    "base_editor_type"
    ]
    if not args.be_files:
        return None
    else:
        ext = Path(args.be_files).suffix.lower()

    if ext not in [".csv", ".tsv", ".txt"]:
        raise ValueError("Unsupported file extension for base editor file. Use .csv, .tsv, or .txt")
    if ext in [".csv"]:
        be_df = pd.read_csv(args.be_files, header=0)
    elif ext in [".tsv", ".txt"]:
        be_df = pd.read_csv(args.be_files, sep="\t", header=0)

    # 列名が期待通りかチェック、違うならエラーを投げる
    if set(be_df.columns) != set(expected_columns):
        raise parser.error(
            f"Base editor file columns are invalid. "
            f"Expected columns: {expected_columns}, but got: {list(be_df.columns)}"
        )
    else:
        return {
            row["base_editor_name"]: BaseEditor(
                base_editor_name=row["base_editor_name"],
                pam_sequence=row["pam_sequence"],
                editing_window_start_in_grna=int(row["editing_window_start"]),
                editing_window_end_in_grna=int(row["editing_window_end"]),
                base_editor_type=row["base_editor_type"],
            )
            for _, row in be_df.iterrows()
        }

def show_base_editors_info(base_editors: dict[str, BaseEditor], parser: argparse.ArgumentParser):
    if base_editors is None:
        parser.error("No base editors available to display.")
    for base_editor in base_editors.values():
        print(f"  - {base_editor.base_editor_name} (Type: {base_editor.base_editor_type}, PAM: {base_editor.pam_sequence}, "
            f"Window: {base_editor.editing_window_start_in_grna}-{base_editor.editing_window_end_in_grna})")
    return None

def parse_base_editors_from_args(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
    base_editors: dict[str, BaseEditor]
) -> dict[str, BaseEditor] | None:
    """
    コマンドライン引数からBaseEditor情報を取得し、base_editorsに追加する。
    必要な情報が揃っていない場合はparser.errorで終了。
    """
    be_name = getattr(args, "be_name", None) or getattr(args, "base_editor_name", None)
    be_pam = getattr(args, "be_pam", None) or getattr(args, "base_editor_pam", None)
    be_window_start = getattr(args, "be_window_start", None) or getattr(args, "base_editor_window_start", None)
    be_window_end = getattr(args, "be_window_end", None) or getattr(args, "base_editor_window_end", None)
    be_type = getattr(args, "be_type", None) or getattr(args, "base_editor_type", None)
    # どれか一つでも指定されていれば、全て必須
    if any([be_name, be_pam, be_window_start, be_window_end, be_type]):
        if not all([be_name, be_pam, be_window_start, be_window_end, be_type]):
            parser.error("Base editor information is incomplete. Please provide all required parameters.")
        base_editors[be_name] = BaseEditor(
            base_editor_name=be_name,
            pam_sequence=be_pam.upper(),
            editing_window_start_in_grna=int(be_window_start),
            editing_window_end_in_grna=int(be_window_end),
            base_editor_type=be_type.lower(),
        )
    return base_editors

def parse_base_editors_from_all_sources(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser
) -> dict[str, BaseEditor] | None:
    """
    コマンドライン引数からBaseEditor情報を取得し、BaseEditorの辞書を返す。
    ファイル、プリセット、個別指定の順で情報を取得する。
    """
    # デフォルトですべてのプリセットを使うことにする
    base_editors = PRESET_BASE_EDITORS.copy()
    # ファイルからの読み込み
    if args.be_files:
        base_editors.update(
            parse_base_editors_from_file(args, parser) or {}
        )
    # 個別指定からの読み込み
    base_editors = parse_base_editors_from_args(args, parser, base_editors) or base_editors
    # 選択されたbase editorの情報を表示
    logging.info("Designing sgRNAs for the following base editors:")
    show_base_editors_info(base_editors, parser)
    return base_editors

def parse_arguments(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser
) -> tuple[Path, Path, Path, list[str], dict[str, BaseEditor], str]:
    """
    このモジュールに含まれるすべての関数のラッパー関数。
    """
    refflat_path, gtf_path, fasta_path, output_directory = parse_path_from_args(args)
    interest_gene_list = parse_genes_from_args(args, parser)
    base_editors = parse_base_editors_from_all_sources(args, parser)
    assembly_name = parse_assembly_name_from_args(args)
    return refflat_path, gtf_path, fasta_path, output_directory, interest_gene_list, base_editors, assembly_name