import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
from . import (
    for_cli_setting,
    refflat_preprocessor,
    sequence_annotator,
    splicing_event_classifier,
    target_exon_extractor,
    sgrna_designer,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
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
    dir_group = parser.add_argument_group("Input/Output Options")
    dir_group.add_argument(
        "--r", "--refflat-path",
        required=True,
        help="Path of refflat file"
    )
    dir_group.add_argument(
        "--f", "--fasta-path",
        required=True,
        help="Path of FASTA file"
    )
    dir_group.add_argument(
        "--o", "--output-directory",
        required=True,
        help="Directory of the output files"
    )
    gene_group = parser.add_argument_group("Gene Options")
    gene_group.add_argument(
        "--g", "--interest-genes",
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

    refflat_path = Path(args.refflat_path)
    fasta_path = Path(args.fasta_path)
    output_directory = Path(args.output_directory)

    try:
        for_cli_setting.check_input_output_directories(refflat_path, fasta_path, output_directory)
    except (NotADirectoryError, FileNotFoundError) as e:
        logging.error(f"Error with input/output directories: {e}")
        sys.exit(1)

    interest_gene_list = args.interest_genes
    if not interest_gene_list:
        logging.error("No interest genes provided.")
        sys.exit(1)

    preset_base_editors = sgrna_designer.make_preset_base_editors()

    # BaseEditorの決定
    base_editors = []
    if args.be_f:
        try:
            base_editors.extend(for_cli_setting.get_base_editors_from_args(args))
        except ValueError as e:
            logging.error(f"Error getting base editors from args: {e} - please check the file format and content.")
            sys.exit(1)

    if args.be_pre and args.be_pre not in preset_base_editors:
        logging.error(f"Unknown base editor preset: {args.be_pre} - available presets are {list(preset_base_editors.keys())}.")
        sys.exit(1)
    else:
        base_editors.extend(preset_base_editors[args.be_pre])

    if args.be_n or args.be_p or args.be_ws or args.be_we or args.be_t:
        base_editors.extend(for_cli_setting.parse_base_editors(args))

    if not base_editors:
        logging.error("No base editor specified. Please provide at least one base editor. Exiting.")
        sys.exit(1)

    logging.info("Designing sgRNAs for the following base editors:")
    for_cli_setting.show_base_editors_info(base_editors)

    logging.info(f"Using this FASTA file as reference genome: {fasta_path}")

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

    print("running processing of refFlat file...")
    refflat = refflat.drop_duplicates(subset=["name"], keep=False)
    try:
        refflat = refflat_preprocessor.preprocess_refflat(refflat, interest_gene_list)
    except ValueError as e:
        logging.error(f"Error preprocessing refFlat file: {e}")
        sys.exit(1)

    logging.info("Classifying splicing events...")
    classified_refflat = splicing_event_classifier.classify_splicing_events(refflat)
    classified_refflat.to_pickle(output_directory / "processed_refFlat.pickle")
    del refflat

    print("Extracting target exons...")
    try:
        target_exon_df, splice_acceptor_single_exon_df, splice_donor_single_exon_df = target_exon_extractor.wrap_extract_target_exon(classified_refflat)
    except ValueError as e:
        logging.error(f"Error extracting target exons: {e}")
        sys.exit(1)

    logging.info("Annotating sequences to dataframe from genome FASTA...")
    target_exon_df_with_acceptor_and_donor_sequence = sequence_annotator.annotate_sequence_to_splice_sites(
        target_exon_df, splice_acceptor_single_exon_df, splice_donor_single_exon_df, fasta_path
    )
    del splice_acceptor_single_exon_df, splice_donor_single_exon_df

    print("designing sgRNAs...")
    target_exon_df_with_sgrna = sgrna_designer.design_sgrna_for_base_editors(
        target_exon_df=target_exon_df_with_acceptor_and_donor_sequence,
        base_editors=base_editors
    )

    target_exon_df_with_sgrna.to_pickle(output_directory / "target_exon_df_with_sgrna.pickle")

if __name__ == "__main__":
    main()
