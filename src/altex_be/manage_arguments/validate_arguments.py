from altex_be.sgrna_designer import BaseEditor
import pandas as pd
from pathlib import Path

def show_base_editors_info(base_editors: dict[str, BaseEditor]):
    if base_editors is None:
        raise ValueError("No base editors available to display.")

    for base_editor in base_editors.values():
        print(f"  - {base_editor.base_editor_name} (Type: {base_editor.base_editor_type}, PAM: {base_editor.pam_sequence}, "
            f"Window: {base_editor.editing_window_start_in_grna}-{base_editor.editing_window_end_in_grna})")

def check_input_output_directories(refflat_path: Path, fasta_path: Path, output_directory: Path):
    if not refflat_path.is_file():
        raise FileNotFoundError(f"The provided refFlat file '{refflat_path}' does not exist.")
    if not fasta_path.is_file():
        raise FileNotFoundError(f"The provided FASTA file '{fasta_path}' does not exist.")
    if not output_directory.is_dir():
        raise NotADirectoryError(f"The provided output directory '{output_directory}' does not exist.")

def load_supported_assemblies() -> list[str]:
    """
    パッケージ内のcrispr_direct_supported_assemblies.txtを読み込み、アセンブリ名リストを返す
    """
    # このファイルと同じディレクトリにあるtxtを参照
    txt_path = Path(__file__).parent / "crispr_direct_supported_assemblies.txt"
    with open(txt_path, encoding="utf-8") as f:
        supported_assemblies = {line.strip() for line in f if line.strip() and not line.startswith("#")}
    return supported_assemblies

def is_supported_assembly_name_in_crispr_direct(assembly_name:str) -> bool:
    """
    Purpose: ユーザーが入力したアセンブリ名がCRISPRdirectでサポートされているかを確認する。
    Parameter: 
            assembly_name (str): ユーザーが入力したアセンブリ名
            supported_assemblies (list[str]): CRISPRdirectでサポートされているアセンブリ名のリスト 長いので外部txtとして保存
    return : bool
    """
    supported_assemblies = load_supported_assemblies()
    if assembly_name not in supported_assemblies:
        return False
    return True

def split_df_by_column_chunks(df: pd.DataFrame, chunk_sizes=[12, 6, 6]) -> list[pd.DataFrame]:
    """
    DataFrameをchunk_sizesで指定したカラム数ごとに分割し、各DataFrameをリストで返す。
    """
    columns = df.columns.tolist()
    idx = 0
    df = df.head()
    dfs = []
    for size in chunk_sizes:
        cols = columns[idx:idx+size]
        if not cols:
            break
        dfs.append(df[cols].copy())
        idx += size
    # 残りのカラムも追加
    if idx < len(columns):
        cols = columns[idx:]
        dfs.append(df[cols].copy())
    return dfs