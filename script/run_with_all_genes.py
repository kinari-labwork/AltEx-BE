import pandas as pd
from pathlib import Path
from altex_be import (
    refflat_preprocessor,
    splicing_event_classifier,
    target_exon_extractor,
    sequence_annotator,
    sgrna_designer,
    output_formatter,
    offtarget_scorer,
)
from altex_be.class_def.base_editors import PRESET_BASE_EDITORS, BaseEditor

assembly_names = [
    "hg38", # human
    "mm39", # mouse
    "rn7",  # rat
    "danRer11" # zebrafish
]

def preprocess_refflat_with_all_genes(refflat: pd.DataFrame, gtf_flag: bool) -> pd.DataFrame:
    """
    メインパイプラインと違って、すべての遺伝子を対象にrefflatを前処理する関数
    """
    if refflat.empty:
        return refflat
    refflat = refflat_preprocessor.parse_exon_coordinates(refflat)
    refflat = refflat_preprocessor.calculate_exon_lengths(refflat)
    refflat = refflat_preprocessor.drop_abnormal_mapped_transcripts(refflat)
    refflat = refflat_preprocessor.annotate_coding_information(refflat, gtf_flag)
    refflat = refflat_preprocessor.annotate_frame_information(refflat)
    refflat = refflat_preprocessor.add_exon_position_flags(refflat)
    refflat = refflat_preprocessor.annotate_utr_and_cds_exons(refflat)
    refflat = refflat_preprocessor.add_common_exon_window(refflat)
    refflat = refflat_preprocessor.flag_outside_common_exon_space(refflat)

    return refflat

def run_with_all_genes(
    assembly_name: str
):
    in_output_dir = f"data/{assembly_name}"
    print(f"Processing assembly: {assembly_name}")

    out_file = Path(in_output_dir) / "sgrna_designed_all_genes.pkl"
    formatted_out_file = Path(in_output_dir) / "formatted_exploded_sgrna_all_genes.pkl"
    offtarget_out_file = Path(in_output_dir) / "exploded_sgrna_with_offtarget_all_genes.pkl"
    
    base_editors = dict(PRESET_BASE_EDITORS)
    base_editors.update({
        "spry_be4max_nnn": BaseEditor(
            base_editor_name="spry_be4max_nnn",
            pam_sequence="NNN",
            editing_window_start_in_grna=12,
            editing_window_end_in_grna=17,
            base_editor_type="cbe"
        )})

    refflat = pd.read_csv(
        in_output_dir + "/refFlat.txt",
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
    refflat_preprocessed = preprocess_refflat_with_all_genes(refflat, gtf_flag=False)
    refflat_preprocessed.to_pickle(f"{in_output_dir}/processed_refflat_all_genes.pkl")
    refflat_classified = splicing_event_classifier.classify_splicing_events(refflat_preprocessed)
    refflat_classified.to_pickle(f"{in_output_dir}/classified_refflat_all_genes.pkl")

    splice_acceptor_single_exon_df, splice_donor_single_exon_df, exploded_classified_refflat = target_exon_extractor.wrap_extract_target_exon(refflat_classified)
    exploded_classified_refflat.to_pickle(f"{in_output_dir}/exploded_classified_refflat_all_genes.pkl")

    splice_acceptor_single_exon_df =  splice_acceptor_single_exon_df[splice_acceptor_single_exon_df["chromStart"] > 0]
    sequence_annotated = sequence_annotator.annotate_sequence_to_splice_sites(
        exploded_classified_refflat,
        splice_acceptor_single_exon_df,
        splice_donor_single_exon_df,
        fasta_path=f"{in_output_dir}/{assembly_name}.fa",
    )

    sgrna_designed = sgrna_designer.design_sgrna_for_base_editors(sequence_annotated, base_editors)
    sgrna_designed.to_pickle(out_file)

    sgrna_designed_dict = sgrna_designer.design_sgrna_for_base_editors_dict(sequence_annotated, base_editors)
    formatted_exploded_sgrna_df = output_formatter.format_output(sgrna_designed_dict, base_editors)
    formatted_exploded_sgrna_df.to_pickle(formatted_out_file)

    exploded_sgrna_with_offtarget_info = offtarget_scorer.score_offtargets(
        formatted_exploded_sgrna_df,
        assembly_name,
        fasta_path=Path(f"{in_output_dir}/{assembly_name}.fa"),
    )
    exploded_sgrna_with_offtarget_info.to_pickle(offtarget_out_file)

if __name__ == "__main__":
    for assembly_name in assembly_names:
        run_with_all_genes(assembly_name)
