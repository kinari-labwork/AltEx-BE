import pandas as pd
from pathlib import Path
from altex_be import refflat_preprocessor,splicing_event_classifier, target_exon_extractor, sequence_annotator, sgrna_designer
from altex_be.class_def.base_editors import PRESET_BASE_EDITORS

assembly_names = [
    "hg38", # human
    "mm39", # mouse
    "rn6",  # rat
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
    refflat = refflat_preprocessor.annotate_cording_information(refflat, gtf_flag)
    refflat = refflat_preprocessor.annotate_flame_information(refflat)
    refflat = refflat_preprocessor.add_exon_position_flags(refflat)
    refflat = refflat_preprocessor.annotate_utr_and_cds_exons(refflat)

    return refflat

def run_with_all_genes(
    assembly_name: str
):
    in_output_dir = f"data/{assembly_name}"
    print(f"Processing assembly: {assembly_name}")

    out_file = Path(in_output_dir) / "sgrna_designed_all_genes.pkl"

    if out_file.exists():
        print(f"Skipping {assembly_name}, output already exists: {out_file}")
        return
    
    base_editors = PRESET_BASE_EDITORS

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

if __name__ == "__main__":
    for assembly_name in assembly_names:
        run_with_all_genes(assembly_name)