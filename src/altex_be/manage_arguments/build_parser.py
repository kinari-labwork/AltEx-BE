import argparse
import importlib.metadata

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Altex BE: A CLI tool for processing refFlat files and extracting target exons.",
    )    
    # 明示的に -v/--version を追加
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=importlib.metadata.version("altex-be"),
        help="Show the version of Altex BE",
    )
    # コマンドライン引数を追加
    dir_group = parser.add_argument_group("Input/Output Options")
    dir_group.add_argument(
        "-r", "--refflat-path",
        required=True,
        help="Path of refflat file"
    )
    dir_group.add_argument(
        "-f", "--fasta-path",
        required=True,
        help="Path of FASTA file"
    )
    dir_group.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Directory of the output files"
    )
    gene_group = parser.add_argument_group("Gene Options")
    gene_group.add_argument(
        "--gene-symbols",
        default=None,
        nargs="+",
        help="List of interest gene symbols (space-separated)"
    )
    gene_group.add_argument(
        "--refseq-ids",
        default=None,
        nargs="+",
        help="List of interest gene Refseq IDs (space-separated)"
    )
    gene_group.add_argument(
        "-a", "--assembly-name",
        default=None,
        required=True,
        help="Name of the genome assembly to use"
    )
    gene_group.add_argument(
        "--gene-file",
        default=None,
        required=False,
        help="Path to a file (csv,txt,tsv) containing gene symbols or Refseq IDs (one gene per line)"
    )
    base_editors = parser.add_argument_group("Base Editor Options")
    base_editors.add_argument(
        "-n", "--be-name",
        default=None,
        required=False,
        help="Name of the base editor to optional use",
    )
    base_editors.add_argument(
        "-p", "--be-pam",
        default=None,
        required=False,
        help="PAM sequence for the base editor",
    )
    base_editors.add_argument(
        "-s", "--be-start",
        default=None,
        required=False,
        help="Window start for the base editor (Count from next to PAM)",
    )
    base_editors.add_argument(
        "-e", "--be-end",
        default=None,
        required=False,
        help="Window end for the base editor (Count from next to PAM)",
    )
    base_editors.add_argument(
        "-t", "--be-type",
        default=None,
        required=False,
        help="Choose the type of base editor, this tool supports ABE and CBE",
    )
    base_editors.add_argument(
        "--be-preset",
        default=None,
        required=False,
        help="Preset for the base editor",
    )
    base_editors.add_argument(
        "--be-files",
        default=None,
        required=False,
        help="input the path of csv file or txt file of base editor information",
    )
    return parser