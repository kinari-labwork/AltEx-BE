import argparse
import pandas as pd
from pathlib import Path
import logging
import datetime
from . import (
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

    refflat = loading_and_preprocess_refflat(refflat_path, interest_gene_list, parser)

    classified_refflat = classify_splicing_events(refflat)
    del refflat

    splice_acceptor_single_exon_df, splice_donor_single_exon_df, exploded_classified_refflat = extract_target_exon(
        classified_refflat, interest_gene_list, parser
    )
    del classified_refflat

    target_exon_df_with_acceptor_and_donor_sequence = annotate_sequences(
        exploded_classified_refflat,
        splice_acceptor_single_exon_df,
        splice_donor_single_exon_df,
        fasta_path
    )
    del splice_acceptor_single_exon_df, splice_donor_single_exon_df

    target_exon_df_with_sgrna_dict = design_sgrnas(
        target_exon_df_with_acceptor_and_donor_sequence,
        base_editors
    )
    formatted_exploded_sgrna_df = format_output(target_exon_df_with_sgrna_dict, base_editors, parser)
    del target_exon_df_with_acceptor_and_donor_sequence, exploded_classified_refflat
    
    exploded_sgrna_with_offtarget_info = score_offtargets(formatted_exploded_sgrna_df, assembly_name, fasta_path=fasta_path)

    output_track_name = f"{datetime.datetime.now().strftime('%Y%m%d%H%M')}_{assembly_name}_sgrnas_designed_by_altex-be"
    save_results(exploded_sgrna_with_offtarget_info, output_directory, output_track_name)
    write_ucsc_custom_track(
        exploded_sgrna_with_offtarget_info,
        output_directory,
        output_track_name,
        assembly_name
    )
    
    logging.info("All AltEx-BE processes completed successfully.")
    logging.info("Printing summary of output:")
    summary_dfs = split_df_by_column_chunks(exploded_sgrna_with_offtarget_info, chunk_sizes=[12, 6, 6])
    for sub_df in summary_dfs:
        print(sub_df)  # indexも表示される
        print("-" * 40)
    return

def loading_and_preprocess_refflat(refflat_path: str, interest_gene_list: list[str], parser: argparse.ArgumentParser) -> pd.DataFrame:
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
    refflat = refflat_preprocessor.preprocess_refflat(refflat, interest_gene_list)
    if refflat.empty :
        parser.error("No interest genes found in refFlat after preprocessing. Exiting...")
    # すべて constitutive exonでも設計対象とするが、exonが1つしかない遺伝子は対象外とする
    if not refflat_preprocessor.check_multiple_exon_existance(refflat, interest_gene_list) :
        parser.error("all of your interest genes are single-exon genes. AltEx-BE cannot process these genes. Exiting...")
    return refflat

def classify_splicing_events(refflat: pd.DataFrame) -> pd.DataFrame:
    """
    前処理されたrefflatデータフレームから、スプライシングイベントを分類する。
    """
    logging.info("-" * 50)
    logging.info("Classifying splicing events...")
    classified_refflat = splicing_event_classifier.classify_splicing_events(refflat)
    return classified_refflat

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

def annotate_sequences(
    exploded_classified_refflat: pd.DataFrame,
    splice_acceptor_single_exon_df: pd.DataFrame,
    splice_donor_single_exon_df: pd.DataFrame,
    fasta_path: str
) -> pd.DataFrame:
    logging.info("-" * 50)
    logging.info("Annotating sequences to dataframe from genome FASTA...")
    logging.info(f"Using this FASTA file as reference genome: {fasta_path}")
    target_exon_df_with_acceptor_and_donor_sequence = sequence_annotator.annotate_sequence_to_splice_sites(
        exploded_classified_refflat, splice_acceptor_single_exon_df, splice_donor_single_exon_df, fasta_path
    )
    return target_exon_df_with_acceptor_and_donor_sequence

def design_sgrnas(
    target_exon_df_with_acceptor_and_donor_sequence: pd.DataFrame,
    base_editors: dict[str, BaseEditor]
) -> dict[str, pd.DataFrame]:
    logging.info("designing sgRNAs...")
    target_exon_df_with_sgrna_dict = sgrna_designer.design_sgrna_for_base_editors_dict(
        target_exon_df=target_exon_df_with_acceptor_and_donor_sequence,
        base_editors=base_editors
    )
    return target_exon_df_with_sgrna_dict

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

def score_offtargets(
    formatted_exploded_sgrna_df: pd.DataFrame,
    assembly_name: str,
    fasta_path: str
) -> pd.DataFrame:
    logging.info("-" * 50)
    logging.info("Scoring off-targets...")
    exploded_sgrna_with_offtarget_info = offtarget_scorer.score_offtargets(formatted_exploded_sgrna_df, assembly_name, fasta_path=fasta_path)
    logging.info("-" * 50)
    return exploded_sgrna_with_offtarget_info

def save_results(
    exploded_sgrna_with_offtarget_info: pd.DataFrame,
    output_directory: Path,
    output_track_name: str
) -> None:
    logging.info("Saving results...")
    exploded_sgrna_with_offtarget_info.to_csv(output_directory / f"{output_track_name}_table.csv")
    logging.info(f"Results saved to: {output_directory / f'{output_track_name}_table.csv'}")
    return

def write_ucsc_custom_track(
    exploded_sgrna_with_offtarget_info: pd.DataFrame,
    output_directory: Path,
    output_track_name: str,
    assembly_name: str
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

if __name__ == "__main__":
    run_pipeline()
