import argparse
import pandas as pd
from pathlib import Path
from . import validate_arguments
from .. class_def.base_editors import BaseEditor, PRESET_BASE_EDITORS

def parse_path_from_args(args: argparse.Namespace):
    refflat_path = Path(args.refflat_path)
    fasta_path = Path(args.fasta_path)
    output_directory = Path(args.output_dir)
    return refflat_path, fasta_path, output_directory

def parse_genes_from_args(args: argparse.Namespace):
    genes_from_file = validate_arguments.parse_gene_file(Path(args.gene_file)) if args.gene_file else []
    gene_symbols = args.gene_symbols if args.gene_symbols is not None else []
    refseq_ids = args.refseq_ids if args.refseq_ids is not None else []
    interest_gene_list = gene_symbols + refseq_ids + genes_from_file
    return interest_gene_list

def parse_base_editors_from_file(
    args: argparse.Namespace, 
    parser: argparse.ArgumentParser, 
    base_editors: dict[str, BaseEditor]
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
    if not args.base_editor_files:
        return None
    else:
        ext = Path(args.base_editor_files).suffix.lower()

    if ext not in [".csv", ".tsv", ".txt"]:
        raise ValueError("Unsupported file extension for base editor file. Use .csv, .tsv, or .txt")
    if ext in [".csv"]:
        be_df = pd.read_csv(args.base_editor_files, header=0)
    elif ext in [".tsv", ".txt"]:
        be_df = pd.read_csv(args.base_editor_files, sep=None, engine="python", header=0)

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

def parse_base_editors_from_presets(
    args: argparse.Namespace, 
    parser: argparse.ArgumentParser, 
    base_editors: dict[str, BaseEditor]
) -> dict[str, BaseEditor] | None:
    preset_base_editors = PRESET_BASE_EDITORS # 事前定義されたBaseEditorの辞書
    if args.be_preset is not None:
        if args.be_preset not in preset_base_editors:
            parser.error(f"Invalid base editor preset: {args.be_preset}. Available presets are: {list(preset_base_editors.keys())}")
        else:
            base_editor = preset_base_editors[args.be_preset]
            base_editors[base_editor.base_editor_name] = base_editor
    return base_editors

def parse_base_editors_from_args(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
    base_editors: dict[str, BaseEditor]
) -> dict[str, BaseEditor] | None:
    """
    コマンドライン引数からBaseEditor情報を取得し、base_editorsに追加する。
    必要な情報が揃っていない場合はparser.errorで終了。
    """
    required_fields = [
        args.base_editor_name,
        args.base_editor_pam,
        args.base_editor_window_start,
        args.base_editor_window_end,
        args.base_editor_type
    ]
    # どれか一つでも指定されていれば、全て必須
    if any([args.base_editor_name, args.base_editor_pam, args.base_editor_window_start, args.base_editor_window_end, args.base_editor_type]):
        if not all(required_fields):
            parser.error("Base editor information is incomplete. Please provide all required parameters.")
        base_editors[args.base_editor_name] = BaseEditor(
            base_editor_name=args.base_editor_name,
            pam_sequence=args.base_editor_pam.upper(),
            editing_window_start_in_grna=int(args.base_editor_window_start),
            editing_window_end_in_grna=int(args.base_editor_window_end),
            base_editor_type=args.base_editor_type.lower(),
        )
    return base_editors