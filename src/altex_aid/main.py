import argparse
import pandas as pd
from pathlib import Path
import sys
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
    dir_group = parser.add_argument_group("Directory Options")
    dir_group.add_argument(
        "-i", "--input-directory",
        required=True,
        help="Directory of the input files"
    )
    dir_group.add_argument(
        "-o", "--output-directory",
        required=True,
        help="Directory of the output files"
    )
    gene_group = parser.add_argument_group("Gene Options")
    gene_group.add_argument(
        "-g", "--interest-genes",
        nargs="+",
        required=True,
        help="List of interest genes (space-separated)"
    )
    base_editors = parser.add_argument_group("Base Editor Options")
    base_editors.add_argument(
        "--be-n", "--base-editor-name",
        default=None,
        required=False,
        help="Name of the base editor to optional use",
    )
    base_editors.add_argument(
        "--be-p", "--base-editor-pam",
        default=None,
        required=False,
        help="PAM sequence for the base editor",
    )
    base_editors.add_argument(
        "--be-ws", "--base-editor-window-start",
        default=None,
        required=False,
        help="Window start for the base editor (Count from next to PAM)",
    )
    base_editors.add_argument(
        "--be-we", "--base-editor-window-end",
        default=None,
        required=False,
        help="Window end for the base editor (Count from next to PAM)",
    )
    base_editors.add_argument(
        "--be-t", "--base-editor-type",
        default=None,
        required=False,
        help="Choose the type of base editor, this tool supports ABE and CBE",
    )
    base_editors.add_argument(
        "--be-pre", "--base-editor-preset",
        default=None,
        required=False,
        help="Preset for the base editor",
    )
    base_editors.add_argument(
        "--be-f", "--base-editor-files",
        default=None,
        required=False,
        help="input the path of csv file or txt file of base editor information",
    )

    args = parser.parse_args()

    input_directory = Path(args.input_directory)
    output_directory = Path(args.output_directory)

    try:
        for_cli_setting.check_input_output_directories(input_directory, output_directory)
    except (NotADirectoryError, FileNotFoundError) as e:
        print(e)
        sys.exit(1)

    interest_gene_list = args.interest_genes
    if not interest_gene_list:
        print("No interest genes provided.")
        sys.exit(1)

    preset_base_editors = sgrna_designer.make_preset_base_editors()

    # BaseEditorの決定
    if args.be_f:
        base_editors = for_cli_setting.get_base_editors_from_args(args)
    elif args.be_pre and args.be_pre in preset_base_editors:
        base_editors = [preset_base_editors[args.be_pre]]
    else:
        base_editors = for_cli_setting.parse_base_editors(args)
        if not base_editors:
            print("No base editor specified.")
            sys.exit(1)

    print("Designing sgRNAs for the following base editors:")
    for_cli_setting.show_base_editors_info(base_editors)

    # FASTA ファイルの検出と確認
    fasta_files = [f for f in input_directory.glob("*.fa")] + [f for f in input_directory.glob("*.fasta")]
    # 1 つだけ存在する場合、そのファイルを使用
    try:
        for_cli_setting.check_fasta_files(input_directory, fasta_files)
    except (FileNotFoundError, FileExistsError) as e:
        print(e)
        sys.exit(1)

    fasta_path = input_directory / fasta_files[0]
    print(f"Using this FASTA file as reference genome: {fasta_path}")

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
    try:
        refflat = refflat_preprocessor.preprocess_refflat(refflat, interest_gene_list)
    except ValueError as e:
        print(e)
        sys.exit(1)

    print("Classifying splicing events...")
    classified_refflat = splicing_event_classifier.classify_splicing_events(refflat)
    classified_refflat.to_pickle(f"{output_directory}/processed_refFlat.pickle")
    del refflat

    print("Extracting target exons...")
    try:
        target_exon_df, splice_acceptor_single_exon_df, splice_donor_single_exon_df = target_exon_extractor.wrap_extract_target_exon(classified_refflat)
    except ValueError as e:
        print(e)
        sys.exit(1)

    print("Annotating sequences to dataframe from genome FASTA...")
    target_exon_df_with_acceptor_and_donor_sequence = sequence_annotator.annotate_sequence_to_splice_sites(
        target_exon_df, splice_acceptor_single_exon_df, splice_donor_single_exon_df, fasta_path
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
