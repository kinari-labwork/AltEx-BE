import pandas as pd
from altex_be.output_formatter import (
    convert_empty_list_into_na,
    prepare_melted_df,
    explode_sgrna_df,
    add_base_editor_info_to_df
)
from altex_be.sgrna_designer import BaseEditor


def test_convert_empty_list_into_na():
    input_df = pd.DataFrame({
        "A": [1, 2, []],
        "B": [[], 5, 6]
    })
    expected_df = pd.DataFrame({
        "A": [1, 2, pd.NA],
        "B": [pd.NA, 5, 6]
    })
    output_dict = convert_empty_list_into_na({"test": input_df})
    pd.testing.assert_frame_equal(output_dict["test"], expected_df)

def test_prepare_melted_df():
    input_df = pd.DataFrame({
        "geneName": ["gene1", "gene2"],
        "chrom": ["chr1", "chr1"],
        "exonStarts": [100, 200],
        "exonEnds": [150, 250],
        "strand": ["+", "-"],
        "exonlengths": [50, 50],
        "exontype": ["alternative", "unique-alternative"],
        "coding": ["coding", "coding"],
        "frame" : ["out-frame", "in-frame"],
        "exon_position": ["first", "last"],
        "uuid" : ["uuid1", "uuid2"],
        "cds_info": ["utr_exon", "cds_exon"],
        "acceptor_sgrna_target_sequence": [["CCC+ATCG", "CCC+TAGC"], ["CCC+CGTA", "CCC+GCAT"]],
        "donor_sgrna_target_sequence": [["CCC+GCTA", "CCC+TACG"], ["CCC+AGCT", "CCC+CTAG"]]
    })
    expected_df = pd.DataFrame({
        "geneName": ["gene1", "gene2", "gene1", "gene2"],
        "chrom": ["chr1", "chr1", "chr1", "chr1"],
        "exonStarts": [100, 200, 100, 200],
        "exonEnds": [150, 250, 150, 250],
        "strand": ["+", "-", "+", "-"],
        "exonlengths": [50, 50, 50, 50],
        "coding": ["coding", "coding", "coding", "coding"],
        "frame": ["out-frame", "in-frame", "out-frame", "in-frame"],
        "exontype": ["alternative", "unique-alternative", "alternative", "unique-alternative"],
        "exon_position": ["first", "last", "first", "last"],
        "cds_info": ["utr_exon", "cds_exon", "utr_exon", "cds_exon"],
        "uuid": ["uuid1", "uuid2", "uuid1", "uuid2"],
        "sgrna_target_sequence": [["CCC+ATCG", "CCC+TAGC"], ["CCC+CGTA", "CCC+GCAT"], ["CCC+GCTA", "CCC+TACG"], ["CCC+AGCT", "CCC+CTAG"]],
        "site_type": ["acceptor", "acceptor", "donor", "donor"]
    })
    output_df = prepare_melted_df(input_df)
    print(output_df)
    pd.testing.assert_frame_equal(
        output_df,
        expected_df
    )

def test_explode_sgrna_df():
    input_df1 = pd.DataFrame({
        "geneName": ["gene1"],
        "chrom": ["chr1"],
        "exonStarts": [100],
        "exonEnds": [150],
        "strand": ["+"],
        "exonlengths": [50],
        "coding": ["coding"],
        "frame": ["out-frame"],
        "exontype": ["alternative"],
        "exon_position": ["first"],
        "cds_info": ["utr_exon"],
        "uuid": ["uuid1"],
        "acceptor_sgrna_target_sequence": [["CCC+ATCG", "CCC+TAGC"]],
        "donor_sgrna_target_sequence": [["CCC+GCTA", "CCC+TACG"]]
    })
    input_df2 = pd.DataFrame({
        "geneName": ["gene2"],
        "chrom": ["chr1"],
        "exonStarts": [200],
        "exonEnds": [250],
        "strand": ["-"],
        "exonlengths": [50],
        "coding": ["non-coding"],
        "frame": ["in-frame"],
        "exontype": ["unique-alternative"],
        "exon_position": ["last"],
        "cds_info": ["cds_exon"],
        "uuid": ["uuid2"],
        "acceptor_sgrna_target_sequence": [["CCC+CGTA", "CCC+GCAT"]],
        "donor_sgrna_target_sequence": [["CCC+AGCT", "CCC+CTAG"]]
    })

    input_dict = {
        "MockBE1": input_df1,
        "MockBE2": input_df2
    }

    expected_df = pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene1", "gene1", "gene2", "gene2", "gene2", "gene2"],
        "chrom": ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
        "exonStarts": [100, 100, 100, 100, 200, 200, 200, 200],
        "exonEnds": [150, 150, 150, 150, 250, 250, 250, 250],
        "strand": ["+", "+", "+", "+", "-", "-", "-", "-"],
        "exonlengths": [50, 50, 50, 50, 50, 50, 50, 50],
        "coding": ["coding", "coding", "coding", "coding", "non-coding", "non-coding", "non-coding", "non-coding"],
        "frame": ["out-frame", "out-frame", "out-frame", "out-frame", "in-frame", "in-frame", "in-frame", "in-frame"],
        "exontype": ["alternative", "alternative", "alternative", "alternative", "unique-alternative", "unique-alternative", "unique-alternative", "unique-alternative"],
        "exon_position": ["first", "first", "first", "first", "last", "last", "last", "last"],
        "cds_info": ["utr_exon", "utr_exon", "utr_exon", "utr_exon", "cds_exon", "cds_exon", "cds_exon", "cds_exon"],
        "uuid": ["uuid1", "uuid1", "uuid1", "uuid1", "uuid2", "uuid2", "uuid2", "uuid2"],
        "sgrna_target_sequence": ["CCC+ATCG", "CCC+TAGC", "CCC+GCTA", "CCC+TACG", "CCC+CGTA", "CCC+GCAT", "CCC+AGCT", "CCC+CTAG"],
        "site_type": ["acceptor", "acceptor", "donor", "donor", "acceptor", "acceptor", "donor", "donor"],
        "base_editor_name": ["MockBE1", "MockBE1", "MockBE1", "MockBE1", "MockBE2", "MockBE2", "MockBE2", "MockBE2"],
    })
    output_df = explode_sgrna_df(input_dict)
    print(output_df)

    pd.testing.assert_frame_equal(
        output_df.sort_values(by=list(output_df.columns)).reset_index(drop=True),
        expected_df.sort_values(by=list(expected_df.columns)).reset_index(drop=True)
    )

def test_add_base_editor_info_to_df():
    input_df = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr1"],
        "exonStarts": [100, 100, 100, 100],
        "exonEnds": [150, 150, 150, 150],
        "strand": ["+", "+", "+", "+"],
        "exontype": ["alternative", "alternative", "alternative", "alternative"],
        "exon_position": ["first", "first", "first", "first"],
        "base_editor_name": ["MockBE1", "MockBE1", "MockBE2", "MockBE2"],
        "sgrna_target_sequence": ["CCC+ATCG", "CCC+TAGC", "CCC+GCTA", "CCC+TACG"],
        "site_type": ["acceptor", "acceptor", "donor", "donor"]
    })

    mock_base_editors = {
        "MockBE1": BaseEditor(
            base_editor_name="MockBE1",
            pam_sequence="NGG",
            editing_window_start_in_grna=10,
            editing_window_end_in_grna=15,
            base_editor_type="cbe"
        ),
        "MockBE2": BaseEditor(
            base_editor_name="MockBE2",
            pam_sequence="NGG",
        editing_window_start_in_grna=15,
        editing_window_end_in_grna=20,
        base_editor_type="abe"
    )}
    expected_df = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr1"],
        "exonStarts": [100, 100, 100, 100],
        "exonEnds": [150, 150, 150, 150],
        "strand": ["+", "+", "+", "+"],
        "exontype": ["alternative", "alternative", "alternative", "alternative"],
        "exon_position": ["first", "first", "first", "first"],
        "base_editor_name": ["MockBE1", "MockBE1", "MockBE2", "MockBE2"],
        "base_editor_pam_sequence": ["NGG", "NGG", "NGG", "NGG"],
        "base_editor_editing_window_start": [10, 10, 15, 15],
        "base_editor_editing_window_end": [15, 15, 20, 20],
        "base_editor_type": ["cbe", "cbe", "abe", "abe"],
        "sgrna_target_sequence": ["CCC+ATCG", "CCC+TAGC", "CCC+GCTA", "CCC+TACG"],
        "site_type": ["acceptor", "acceptor", "donor", "donor"],
    })
    output_df = add_base_editor_info_to_df(input_df, mock_base_editors)
    output_df = output_df[expected_df.columns]  
    print(output_df)
    pd.testing.assert_frame_equal(output_df, expected_df)