import argparse
import pandas as pd
import sys
import os

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

    input_directory = args.input_directory
    if not os.path.isdir(input_directory):
        print(f"The provided input directory '{input_directory}' does not exist. Exiting.")
        sys.exit(0)

    output_directory = args.output_directory
    if not os.path.isdir(output_directory):
        print(f"The provided output directory '{output_directory}' does not exist. Exiting.")
        sys.exit(0)

    interest_gene_list = args.interest_genes.split()
    if not interest_gene_list:
        print("No interest genes provided. Exiting.")
        sys.exit(0)

    base_editors = sgrna_designer.make_default_base_editors()

    #　コマンドライン引数で BaseEditor 情報が提供されている場合は、解析し、base_editors に追加
    if args.base_editors:
        base_editors = for_cli_setting.parse_base_editors(args.base_editors, base_editors)
    
    print("Designing sgRNAs for the following base editors:")
    for base_editor in base_editors:
        print(f"  - {base_editor.base_editor_name} (Type: {base_editor.base_editor_type}, PAM: {base_editor.pam_sequence}, "
            f"Window: {base_editor.editing_window_start_in_grna}-{base_editor.editing_window_end_in_grna})")

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
    if not os.path.isfile(f"{input_directory}/refFlat.txt"):
        print(f"refFlat file not found at '{input_directory}/refFlat.txt'. Exiting.")
        sys.exit(0)

    print("running processing of refFlat file...")
    refflat = refflat.drop_duplicates(subset=["name"], keep=False)
    refflat = refflat_preprocessor.select_interest_genes(refflat, interest_gene_list)
    if refflat.empty:
        print("No interest genes found in the refFlat file. Exiting.")
        sys.exit(0)

    variant_check = refflat_preprocessor.check_transcript_variant(refflat, interest_gene_list)
    if not variant_check:
        print("No transcript variants found for your interest genes. Exiting.")
        sys.exit(0)
    
    refflat = refflat_preprocessor.parse_exon_coordinates(refflat)
    refflat = refflat_preprocessor.calculate_exon_lengths(refflat)
    refflat = refflat_preprocessor.drop_abnormal_mapped_transcripts(refflat)
    refflat = refflat_preprocessor.annotate_cording_information(refflat)
    refflat = refflat_preprocessor.annotate_flame_information(refflat)
    refflat = refflat_preprocessor.annotate_flame_information(refflat)
    refflat = refflat_preprocessor.add_exon_position_flags(refflat)

    print("Classifying splicing events...")
    classified_refflat = splicing_event_classifier.classify_splicing_events_per_gene(refflat)
    classified_refflat = splicing_event_classifier.flip_a3ss_a5ss_on_minus_strand(classified_refflat)
    classified_refflat.to_pickle(f"{output_directory}/processed_refFlat.pickle")
    del refflat, variant_check

    print("Extracting target exons...")
    target_exon_df = target_exon_extractor.extract_target_exon(classified_refflat)
    if not target_exon_extractor.check_target_exon_existence(target_exon_df):
        print("No target exons found. Exiting.")
        sys.exit(0)
    splice_acceptor_single_exon_df = (
        target_exon_extractor.extract_splice_acceptor_regions(target_exon_df, 25)
    )
    splice_donor_single_exon_df = (
        target_exon_extractor.extract_splice_donor_regions(target_exon_df, 25)
    )

    print("Annotating sequences to dataframe from genome FASTA...")
    fasta_path = f"{input_directory}/combined.fa"

    if not os.path.isfile(fasta_path):
        print(f"FASTA file not found at '{fasta_path}'. Exiting.")
        sys.exit(0)

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
