import argparse
import pandas as pd
from pathlib import Path
import logging
import sys
import datetime
from . import (
    cli_setting,
    refflat_preprocessor,
    sequence_annotator,
    splicing_event_classifier,
    target_exon_extractor,
    sgrna_designer,
    output_formatter,
    offtarget_scorer,
    bed_for_ucsc_custom_track_maker,
    logging_config # noqa: F401
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
        nargs="+",
        help="List of interest gene symbols (space-separated)"
    )
    gene_group.add_argument(
        "--refseq-ids",
        nargs="+",
        help="List of interest gene Refseq IDs (space-separated)"
    )
    gene_group.add_argument(
        "-a", "--assembly-name",
        default=None,
        required=True,
        help="Name of the genome assembly to use"
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

    args = parser.parse_args()

    refflat_path = Path(args.refflat_path)
    fasta_path = Path(args.fasta_path)
    output_directory = Path(args.output_dir)

    cli_setting.check_input_output_directories(refflat_path, fasta_path, output_directory)

    gene_symbols = args.gene_symbols if args.gene_symbols is not None else []
    refseq_ids = args.refseq_ids if args.refseq_ids is not None else []
    interest_gene_list = gene_symbols + refseq_ids
    if not interest_gene_list:
        raise ValueError("Please provide at least one interest gene symbol or Refseq ID.")

    preset_base_editors = sgrna_designer.make_preset_base_editors()

    # BaseEditorの決定
    base_editors = {}
    if args.be_files:
        base_editors = cli_setting.get_base_editors_from_args(args)

    if args.be_preset and args.be_preset not in preset_base_editors:
        raise ValueError(f"Invalid base editor preset: {args.be_preset}. Available presets are: {list(preset_base_editors.keys())}")
    else:
        base_editor = preset_base_editors[args.be_preset]
        base_editors[base_editor.base_editor_name] = base_editor

    if args.be_name or args.be_pam or args.be_start or args.be_end or args.be_type:
        base_editors.update(cli_setting.parse_base_editors(args))

    assembly_name = str(args.assembly_name)
    if not cli_setting.is_supported_assembly_name_in_crispr_direct(assembly_name):
        logging.warning(f"your_assembly : {assembly_name} is not supported by CRISPRdirect. please see <https://crispr.dbcls.jp/doc/>")
    
    output_track_name = f"{datetime.datetime.now().strftime('%Y%m%d%H%M')}_{assembly_name}_sgrnas_designed_by_altex-be"

    if not base_editors:
        raise ValueError("No base editors specified. Please provide at least one base editor.")

    logging.info("Designing sgRNAs for the following base editors:")
    cli_setting.show_base_editors_info(base_editors)

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

    logging.info("running processing of refFlat file...")
    refflat = refflat.drop_duplicates(subset=["name"], keep=False)
    refflat = refflat_preprocessor.preprocess_refflat(refflat, interest_gene_list)
    if not refflat_preprocessor.validate_filtered_refflat(refflat, interest_gene_list) :
        logging.warning("Filtered refFlat DataFrame is invalid. Exiting...")
        sys.exit(0)

    logging.info("Classifying splicing events...")
    classified_refflat = splicing_event_classifier.classify_splicing_events(refflat)
    del refflat

    logging.info("Extracting target exons...")
    splice_acceptor_single_exon_df, splice_donor_single_exon_df, exploded_classified_refflat = target_exon_extractor.wrap_extract_target_exon(classified_refflat)
    if splice_acceptor_single_exon_df.empty and splice_donor_single_exon_df.empty:
        logging.warning("No target exons found for all of the given genes, exiting")
        sys.exit(0)
    for gene in interest_gene_list:
        if gene not in exploded_classified_refflat['geneName'].values:
            logging.info(f"No target exons found for the gene: {gene}. Further processing of {gene} will be skipped.")
        else:
            logging.info(f"Target exons found for the gene: {gene}.")

    logging.info("Annotating sequences to dataframe from genome FASTA...")
    target_exon_df_with_acceptor_and_donor_sequence = sequence_annotator.annotate_sequence_to_splice_sites(
        exploded_classified_refflat, splice_acceptor_single_exon_df, splice_donor_single_exon_df, fasta_path
    )
    del splice_acceptor_single_exon_df, splice_donor_single_exon_df

    logging.info("designing sgRNAs...")
    target_exon_df_with_sgrna_dict = sgrna_designer.design_sgrna_for_base_editors_dict(
        target_exon_df=target_exon_df_with_acceptor_and_donor_sequence,
        base_editors=base_editors
    )

    logging.info("Formatting output...")
    formatted_exploded_sgrna_df = output_formatter.format_output(target_exon_df_with_sgrna_dict, base_editors)
    if formatted_exploded_sgrna_df.empty:
        logging.warning("No sgRNAs could be designed for given genes and Base Editors, Exiting")
        sys.exit(0)
    del target_exon_df_with_acceptor_and_donor_sequence, exploded_classified_refflat

    logging.info("Scoring off-targets...")
    exploded_sgrna_with_offtarget_info = offtarget_scorer.score_offtargets(formatted_exploded_sgrna_df, assembly_name, fasta_path=fasta_path)

    logging.info("Saving results...")
    exploded_sgrna_with_offtarget_info.to_csv(output_directory / f"{output_track_name}_table.csv")

    logging.info("Generating UCSC custom track...")
    bed_df = bed_for_ucsc_custom_track_maker.format_sgrna_for_ucsc_custom_track(exploded_sgrna_with_offtarget_info)

    output_path = output_directory / f"{output_track_name}_ucsc_custom_track.bed"
    track_description: str = f"sgRNAs designed by Altex-BE on {datetime.datetime.now().strftime('%Y%m%d')}"

    with open(output_path, "w") as f:
        track_header = f'track name="{output_track_name}" description="{track_description}" visibility=2 itemRgb="On"\n'
        f.write(track_header)
        bed_df.to_csv(f, sep="\t", header=False, index=False, lineterminator='\n')

    logging.info(f"UCSC custom track file saved to: {output_path}")
    
    return logging.info("All altex-be processes completed successfully.")

if __name__ == "__main__":
    main()
