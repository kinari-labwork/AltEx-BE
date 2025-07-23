import argparse

import pandas as pd

from . import (
    refflat_preprocessor,
    sequence_annotator,
    splicing_event_classifier,
    target_exon_extractor,
    sgrna_designer,
)


def main():
    parser = argparse.ArgumentParser(
        description="Altex AID: A CLI tool for processing refFlat files and extracting target exons.",
    )

    subparsers = parser.add_subparsers(
        dest="mode", required=True, help="Choose a mode: 'run' or 'view'"
    )
    # "run" サブコマンドの定義
    run_parser = subparsers.add_parser("run", help="Run the Altex-AID pipeline.")
    run_parser.add_argument(
        "--input directory",
        required=True,
        help="please input the directory of the input files"
        "This directory should contain the refFlat file and genome FASTA file of your interest species",
    )
    run_parser.add_argument(
        "--output directory",
        required=True,
        help="please input the directory of the output files"
        "This directory will contain the processed refFlat file and other output files",
    )

    # "view" サブコマンドの定義
    view_parser = subparsers.add_parser(
        "view", help="view and search sgRNAs for your interest gene."
    )
    view_parser.add_argument(
        "--output_directory",
        required=True,
        help="please input the directory of the output files"
        "This directory will contain the processed refFlat file and other output files",
    )

    # 明示的に -v/--version を追加
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="0.1.0",
        help="バージョン情報を表示します",
    )

    args = parser.parse_args()

    if args.mode == "run":
        input_directory = args.input_directory
        output_directory = args.output_directory

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
        del classified_refflat

        print("Extracting target exons...")
        target_exon_df = target_exon_extractor.extract_target_exon(refflat)
        splice_acceptor_single_exon_df = (
            target_exon_extractor.extract_splice_acceptor_regions(target_exon_df, 25)
        )
        splice_donor_single_exon_df = (
            target_exon_extractor.extract_splice_donor_regions(target_exon_df, 25)
        )

        print("Annotating sequences to dataframe from genome FASTA...")
        fasta_path = f"{input_directory}/combined.fa"
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
        target_exon_df_with_sgrna_info = sgrna_designer.design_sgrna_for_target_exon_df(
            target_exon_df=target_exon_df_with_acceptor_and_donor_sequence,
            pam_sequence="NGG",
            editing_window_start_in_grna=17,
            editing_window_end_in_grna=19,
        )
        target_exon_df_with_sgrna_info = sgrna_designer.organize_target_exon_df_with_grna_sequence(
            target_exon_df_with_sgrna_info
        )
        target_exon_df_with_sgrna_info = sgrna_designer.convert_sgrna_start_end_position_to_position_in_chromosome(
            target_exon_df_with_sgrna_info
        )


if __name__ == "__main__":
    main()
