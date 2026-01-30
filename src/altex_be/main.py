import argparse
import pandas as pd
from pathlib import Path
import logging
import datetime
import subprocess
import sys
from . import (
    gtf2refflat_converter,
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
from .manage_arguments import (
    build_parser,
    parse_arguments,
    validate_arguments
)
from .class_def.base_editors import BaseEditor


def run_pipeline():
    # --ui flag check before parsing other arguments
    if '--ui' in sys.argv:
        streamlit_app_path = Path(__file__).parent / "streamlit.py"
        if not streamlit_app_path.exists():
            logging.error(f"Streamlit UI file not found at: {streamlit_app_path}")
            sys.exit(1)
        
        logging.info("Launching Streamlit UI...")
        try:
            # Using sys.executable to ensure the same python environment is used
            # Added --server.headless=false to force browser to open
            command = [
                sys.executable, "-m", "streamlit", "run", str(streamlit_app_path),
                "--server.headless=false"
            ]
            subprocess.run(command, check=True)
        except FileNotFoundError:
            logging.error("`streamlit` command not found. Please ensure Streamlit is installed in the current environment.")
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            logging.error(f"Streamlit UI failed to launch: {e}")
            sys.exit(1)
        return # Exit after UI closes

    # --- CLI pipeline ---
    parser = build_parser.build_parser()
    args = parser.parse_args()
    
    refflat_path, gtf_path, fasta_path, output_directory, interest_gene_list, base_editors, assembly_name = parse_arguments.parse_arguments(args, parser)

    validate_arguments.validate_arguments(
        refflat_path,
        gtf_path,
        fasta_path,
        output_directory,
        interest_gene_list,
        base_editors,
        assembly_name,
        parser
    )
    
    if gtf_path is not None :
        logging.info("-" * 50)
        logging.info("Converting GTF to refFlat format...")
        gtf2refflat_converter.gtf_to_refflat(gtf_path, output_directory, assembly_name)
        refflat = loading_and_preprocess_refflat(output_directory / f"converted_refflat_{assembly_name}.txt", interest_gene_list, parser, gtf_flag=True)
    elif refflat_path is not None :
        refflat = loading_and_preprocess_refflat(refflat_path, interest_gene_list, parser, gtf_flag=False)

    logging.info("-" * 50)
    logging.info("Classifying splicing events...")
    classified_refflat = splicing_event_classifier.classify_splicing_events(refflat)
    del refflat

    splice_acceptor_single_exon_df, splice_donor_single_exon_df, exploded_classified_refflat = extract_target_exon(
        classified_refflat, interest_gene_list, parser
    )
    del classified_refflat

    logging.info("-" * 50)
    logging.info("Annotating sequences to dataframe from genome FASTA...")
    logging.info(f"Using this FASTA file as reference genome: {fasta_path}")
    target_exon_df_with_acceptor_and_donor_sequence = sequence_annotator.annotate_sequence_to_splice_sites(
        exploded_classified_refflat, splice_acceptor_single_exon_df, splice_donor_single_exon_df, fasta_path
    )
    del splice_acceptor_single_exon_df, splice_donor_single_exon_df

    logging.info("designing sgRNAs...")
    target_exon_df_with_sgrna_dict = sgrna_designer.design_sgrna_for_base_editors_dict(
        target_exon_df=target_exon_df_with_acceptor_and_donor_sequence,
        base_editors=base_editors
    )

    formatted_exploded_sgrna_df = format_output(target_exon_df_with_sgrna_dict, base_editors, parser)
    del target_exon_df_with_acceptor_and_donor_sequence, exploded_classified_refflat
    
    logging.info("-" * 50)
    logging.info("Scoring off-targets...")
    exploded_sgrna_with_offtarget_info = offtarget_scorer.score_offtargets(formatted_exploded_sgrna_df, assembly_name, fasta_path=fasta_path)
    logging.info("-" * 50)

    output_track_name = f"{datetime.datetime.now().strftime('%Y%m%d%H%M')}_{assembly_name}_sgrnas_designed_by_altex-be"
    logging.info("Saving results...")
    exploded_sgrna_with_offtarget_info.to_csv(output_directory / f"{output_track_name}_table.csv")
    logging.info(f"Results saved to: {output_directory / f'{output_track_name}_table.csv'}")

    write_ucsc_custom_track(
        exploded_sgrna_with_offtarget_info,
        output_directory,
        output_track_name
    )
    
    logging.info("All AltEx-BE processes completed successfully.")
    logging.info(f"Output directory: {output_directory}")
    return

def loading_and_preprocess_refflat(refflat_path: str, interest_gene_list: list[str], parser: argparse.ArgumentParser, gtf_flag: bool) -> pd.DataFrame:
    """
    データのロード、前処理から、興味のある遺伝子の抽出までを行う。
    """
    logging.info("-" * 50)
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
    refflat = refflat_preprocessor.preprocess_refflat(refflat, interest_gene_list, gtf_flag)
    if refflat.empty :
        parser.error("No interest genes found in refFlat after preprocessing. Exiting...")
    # すべて constitutive exonでも設計対象とするが、exonが1つしかない遺伝子は対象外とする
    if not refflat_preprocessor.check_multiple_exon_existance(refflat, interest_gene_list) :
        parser.error("all of your interest genes are single-exon genes. AltEx-BE cannot process these genes. Exiting...")
    return refflat

def extract_target_exon(classified_refflat: pd.DataFrame, interest_gene_list: list[str], parser: argparse.ArgumentParser) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    分類されたスプライシングイベントデータフレームから、ターゲットエキソンを抽出する。
    """
    logging.info("-" * 50)
    logging.info("Extracting target exons...")
    splice_acceptor_single_exon_df, splice_donor_single_exon_df, exploded_classified_refflat = target_exon_extractor.wrap_extract_target_exon(classified_refflat)
    if splice_acceptor_single_exon_df.empty and splice_donor_single_exon_df.empty:
        parser.error("No target exons found for all of the given genes, exiting")
        
    for gene in interest_gene_list:
        if gene not in exploded_classified_refflat['geneName'].values:
            logging.info(f"No target exons found for the gene: {gene}. Further processing of {gene} will be skipped.")
        else:
            logging.info(f"Target exons found for the gene: {gene}.")
    return splice_acceptor_single_exon_df, splice_donor_single_exon_df, exploded_classified_refflat

def format_output(
    target_exon_df_with_sgrna_dict: dict[str, pd.DataFrame],
    base_editors: dict[str, BaseEditor],
    parser: argparse.ArgumentParser
) -> pd.DataFrame:
    logging.info("-" * 50)
    logging.info("Formatting output...")
    formatted_exploded_sgrna_df = output_formatter.format_output(target_exon_df_with_sgrna_dict, base_editors)
    if formatted_exploded_sgrna_df.empty:
        parser.error("No sgRNAs could be designed for given genes and Base Editors, Exiting")
    return formatted_exploded_sgrna_df

def write_ucsc_custom_track(
    exploded_sgrna_with_offtarget_info: pd.DataFrame,
    output_directory: Path,
    output_track_name: str,
) -> None:
    logging.info("Generating UCSC custom track...")
    bed_df = bed_for_ucsc_custom_track_maker.format_sgrna_for_ucsc_custom_track(exploded_sgrna_with_offtarget_info)

    output_path = output_directory / f"{output_track_name}_ucsc_custom_track.bed"
    track_description: str = f"sgRNAs designed by AltEx-BE on {datetime.datetime.now().strftime('%Y%m%d')}"

    with open(output_path, "w") as f:
        track_header = f'track name="{output_track_name}" description="{track_description}" visibility=2 itemRgb="On"\n'
        f.write(track_header)
        bed_df.to_csv(f, sep="\t", header=False, index=False, lineterminator='\n')

    logging.info(f"UCSC custom track file saved to: {output_path}")
    return

if __name__ == "__main__":
    run_pipeline()