import argparse
import pandas as pd
import sys
from pathlib import Path
from . import (
    for_cli_setting,
    refflat_preprocessor,
    sequence_annotator,
    splicing_event_classifier,
    target_exon_extractor,
    sgrna_designer,
)


def main():
    parser = argparse.ArgumentParser(
        description="Altex BE: A CLI tool for processing refFlat files and extracting target exons.",
    )    
    # 明示的に -v/--version を追加
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="0.1.0",
        help="Show the version of Altex BE",
    )
    # コマンドライン引数を追加
    parser.add_argument(
        "-i", "--input-directory",
        required=True,
        help="Directory of the input files"
    )
    parser.add_argument(
        "-o", "--output-directory",
        required=True,
        help="Directory of the output files"
    )
    parser.add_argument(
        "-g", "--interest-genes",
        required=True,
        help="List of interest genes (space-separated)"
    )
    parser.add_argument(
        "-be", "--base-editors",
        nargs='*',
        help=(
            "BaseEditor information in the format: "
            "'name,pam,window_start,window_end,base_editor_type'. "
            "Caution: window_start and window_end are 1-based indices, starting from next to the PAM. "
            "Example: 'ABE8e,NGG,12,17,ABE'"
        ),
    )

    args = parser.parse_args()

    input_directory = Path(args.input_directory)
    if not input_directory.is_dir():
        raise NotADirectoryError(f"The provided input directory '{input_directory}' does not exist.")

    output_directory = Path(args.output_directory)
    if not output_directory.is_dir():
        raise NotADirectoryError(f"The provided output directory '{output_directory}' does not exist.")

    interest_gene_list = args.interest_genes
    if not interest_gene_list:
        raise ValueError("No interest genes provided.")

    base_editors = sgrna_designer.make_default_base_editors()

    #　コマンドライン引数で BaseEditor 情報が提供されている場合は、解析し、base_editors に追加
    if args.base_editors:
        base_editors = for_cli_setting.parse_base_editors(args.base_editors, base_editors)
    
    print("Designing sgRNAs for the following base editors:")

    for base_editor in base_editors:
        print(f"  - {base_editor.base_editor_name} (Type: {base_editor.base_editor_type}, PAM: {base_editor.pam_sequence}, "
            f"Window: {base_editor.editing_window_start_in_grna}-{base_editor.editing_window_end_in_grna})")

    if not (input_directory / "refFlat.txt").is_file():
        raise FileNotFoundError(f"refFlat file not found at '{input_directory}/refFlat.txt'.")

    # FASTA ファイルの検出と確認
    fasta_files = [f for f in input_directory.glob("*.fa")] + [f for f in input_directory.glob("*.fasta")]

    if len(fasta_files) == 0:
        raise FileNotFoundError(f"No FASTA file found in '{input_directory}'.")
    elif len(fasta_files) > 1:
        raise FileExistsError(f"Multiple FASTA files found in '{input_directory}': {', '.join(fasta_files)}. Exiting.")
    # 1 つだけ存在する場合、そのファイルを使用
    fasta_path = input_directory / fasta_files[0]
    print(f"Using FASTA file: {fasta_path}")

    print("loading refFlat file...")
    refflat = pd.read_csv(
            f"{input_directory}/refFlat.txt",
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

    print("running processing of refFlat file...")
    refflat = refflat.drop_duplicates(subset=["name"], keep=False)
    refflat = refflat_preprocessor.preprocess_refflat(refflat, interest_gene_list)

    print("Classifying splicing events...")
    classified_refflat = splicing_event_classifier.classify_splicing_events(refflat)
    classified_refflat.to_pickle(f"{output_directory}/processed_refFlat.pickle")
    del refflat

    print("Extracting target exons...")
    target_exon_df = target_exon_extractor.extract_target_exon(classified_refflat)
    if not target_exon_extractor.check_target_exon_existence(target_exon_df):
        print("No target exons found. Exiting.")
        sys.exit(1)
    splice_acceptor_single_exon_df = (
        target_exon_extractor.extract_splice_acceptor_regions(target_exon_df, 25)
    )
    splice_donor_single_exon_df = (
        target_exon_extractor.extract_splice_donor_regions(target_exon_df, 25)
    )

    print("Annotating sequences to dataframe from genome FASTA...")

    splice_acceptor_single_exon_df = sequence_annotator.annotate_sequence_to_bed(
    splice_acceptor_single_exon_df, fasta_path
    )
    splice_donor_single_exon_df = sequence_annotator.annotate_sequence_to_bed(
        splice_donor_single_exon_df, fasta_path
    )
    target_exon_df_with_acceptor_and_donor_sequence = sequence_annotator.join_sequence_to_single_exon_df(
        single_exon_df=target_exon_df,
        acceptor_bed_with_sequences=splice_acceptor_single_exon_df,
        donor_bed_with_sequences=splice_donor_single_exon_df,
    )
    del splice_acceptor_single_exon_df, splice_donor_single_exon_df

    print("designing sgRNAs...")
    target_exon_df_with_sgrna = sgrna_designer.design_sgrna_for_base_editors(
        target_exon_df=target_exon_df_with_acceptor_and_donor_sequence,
        base_editors=base_editors
    )
    
    target_exon_df_with_sgrna.to_pickle(f"{output_directory}/target_exon_df_with_sgrna.pickle")

if __name__ == "__main__":
    main()
