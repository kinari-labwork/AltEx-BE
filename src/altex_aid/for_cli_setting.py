from altex_aid.sgrna_designer import BaseEditor
import argparse
import pandas as pd
from pathlib import Path


def parse_base_editors(args: argparse.Namespace) -> list[BaseEditor] | None:
    if not all([args.be_n, args.be_p, args.be_ws, args.be_we, args.be_t]):
        raise ValueError(
            "Base editor information is incomplete. Please provide all required parameters."
        )
    try:
        return [
            BaseEditor(
                base_editor_name=args.be_n,
                pam_sequence=args.be_p.upper(),
                editing_window_start_in_grna=int(args.be_ws),
                editing_window_end_in_grna=int(args.be_we),
                base_editor_type=args.be_t.lower(),
            )
        ]
    except ValueError as e:
        print(f"Error parsing base editor information: {e}")
        return None

def show_base_editors_info(base_editors: list[BaseEditor]):
    if base_editors is None:
        print("No base editors available to display.")
        return

    for base_editor in base_editors:
        print(f"  - {base_editor.base_editor_name} (Type: {base_editor.base_editor_type}, PAM: {base_editor.pam_sequence}, "
            f"Window: {base_editor.editing_window_start_in_grna}-{base_editor.editing_window_end_in_grna})")
        

def get_base_editors_from_args(args: argparse.Namespace) -> list[BaseEditor] | None:
    """
    base editorの情報を含むファイルのパスを示す引数を受け取り、BaseEditorのリストを返す。
    csvまたはtxt, tsv形式のファイルをサポートする
    """
    expected_columns = [
    "base_editor_name",
    "pam_sequence",
    "editing_window_start_in_grna",
    "editing_window_end_in_grna",
    "base_editor_type"
    ]
    if not args.be_f:
        return None
    else:
        ext = Path(args.be_f).suffix.lower()

    if ext not in [".csv", ".tsv", ".txt"]:
        raise ValueError("Unsupported file extension for base editor file. Use .csv, .tsv, or .txt")
    if ext in [".csv"]:
        be_df = pd.read_csv(args.be_f, header=0)
    elif ext in [".tsv", ".txt"]:
        be_df = pd.read_csv(args.be_f, sep=None, engine="python", header=0)

    # 列名が期待通りかチェック、違うならエラーを投げる
    if list(be_df.columns) != expected_columns:
        raise ValueError(
            f"Base editor file columns are invalid. "
            f"Expected columns: {expected_columns}, but got: {list(be_df.columns)}"
        )
    else:
        return [
            BaseEditor(
                base_editor_name=row["base_editor_name"],
                pam_sequence=row["pam_sequence"],
                editing_window_start_in_grna=int(row["editing_window_start_in_grna"]),
                editing_window_end_in_grna=int(row["editing_window_end_in_grna"]),
                base_editor_type=row["base_editor_type"],
            )
            for _, row in be_df.iterrows()
        ]

def check_input_output_directories(refflat_path: Path, fasta_path: Path, output_directory: Path):
    if not refflat_path.is_file():
        raise FileNotFoundError(f"The provided refFlat file '{refflat_path}' does not exist.")
    if not fasta_path.is_file():
        raise FileNotFoundError(f"The provided FASTA file '{fasta_path}' does not exist.")
    if not output_directory.is_dir():
        raise NotADirectoryError(f"The provided output directory '{output_directory}' does not exist.")
