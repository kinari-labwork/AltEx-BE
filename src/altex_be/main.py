import argparse
import pandas as pd
from pathlib import Path
import logging
import sys
import datetime
from . import (
    sequence_annotator,
    splicing_event_classifier,
    target_exon_extractor,
    sgrna_designer,
    output_formatter,
    offtarget_scorer,
    bed_for_ucsc_custom_track_maker,
    logging_config # noqa: F401
)
from .manage_arguments import (
    build_parser,
    parse_arguments,
    validate_arguments
)
from .loding_and_preprocess.loding_and_preprocessing import loding_and_preprocess_refflat


def main():
    parser = build_parser.build_parser()

    args = parser.parse_args()

    refflat_path, fasta_path, output_directory, interest_gene_list, base_editors, assembly_name = parse_arguments.parse_arguments(args)

    validate_arguments.validate_arguments(
        refflat_path,
        fasta_path,
        output_directory,
        interest_gene_list,
        base_editors,
        assembly_name,
        parser
    )

    output_track_name = f"{datetime.datetime.now().strftime('%Y%m%d%H%M')}_{assembly_name}_sgrnas_designed_by_altex-be"
    logging.info(f"Using this FASTA file as reference genome: {fasta_path}")

    refflat = loding_and_preprocess_refflat(refflat_path, interest_gene_list, parser)

    logging.info("-" * 50)

    logging.info("Classifying splicing events...")

    classified_refflat = splicing_event_classifier.classify_splicing_events(refflat)
    del refflat
    logging.info("-" * 50)
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
    logging.info("-" * 50)

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
    logging.info("-" * 50)
    logging.info("Formatting output...")
    formatted_exploded_sgrna_df = output_formatter.format_output(target_exon_df_with_sgrna_dict, base_editors)
    if formatted_exploded_sgrna_df.empty:
        logging.warning("No sgRNAs could be designed for given genes and Base Editors, Exiting")
        sys.exit(0)
    del target_exon_df_with_acceptor_and_donor_sequence, exploded_classified_refflat
    
    logging.info("Scoring off-targets...")
    exploded_sgrna_with_offtarget_info = offtarget_scorer.score_offtargets(formatted_exploded_sgrna_df, assembly_name, fasta_path=fasta_path)
    logging.info("-" * 50)
    logging.info("Saving results...")
    exploded_sgrna_with_offtarget_info.to_csv(output_directory / f"{output_track_name}_table.csv")
    logging.info(f"Results saved to: {output_directory / f'{output_track_name}_table.csv'}")

    logging.info("Generating UCSC custom track...")
    bed_df = bed_for_ucsc_custom_track_maker.format_sgrna_for_ucsc_custom_track(exploded_sgrna_with_offtarget_info)

    output_path = output_directory / f"{output_track_name}_ucsc_custom_track.bed"
    track_description: str = f"sgRNAs designed by AltEx-BE on {datetime.datetime.now().strftime('%Y%m%d')}"

    with open(output_path, "w") as f:
        track_header = f'track name="{output_track_name}" description="{track_description}" visibility=2 itemRgb="On"\n'
        f.write(track_header)
        bed_df.to_csv(f, sep="\t", header=False, index=False, lineterminator='\n')

    logging.info(f"UCSC custom track file saved to: {output_path}")
    
    logging.info("All AltEx-BE processes completed successfully.")
    logging.info("Printing summary of output:")
    summary_dfs = validate_arguments.split_df_by_column_chunks(exploded_sgrna_with_offtarget_info, chunk_sizes=[12, 6, 6])
    for sub_df in summary_dfs:
        print(sub_df)  # indexも表示される
        print("-" * 40)
    return

if __name__ == "__main__":
    main()
