import pandas as pd
from altex_aid.offtarget_scorer import (
    convert_empty_list_into_na,
    prepare_melted_df,
    explode_sgrna_df,
    add_base_editor_info_to_df,
    add_crisprdirect_url_to_df,
    calculate_offtarget_site_count_to_df
)
from altex_aid.sgrna_designer import BaseEditor

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
        "chrom": ["chr1", "chr1"],
        "exonStarts": [100, 200],
        "exonEnds": [150, 250],
        "strand": ["+", "-"],
        "exontype": ["skipped", "unique"],
        "exon_position": ["first", "last"],
        "base_editor_name": ["MockBE1", "MockBE2"],
        "acceptor_sgrna_target_sequence": [["CCC+ATCG", "CCC+TAGC"], ["CCC+CGTA", "CCC+GCAT"]],
        "donor_sgrna_target_sequence": [["CCC+GCTA", "CCC+TACG"], ["CCC+AGCT", "CCC+CTAG"]]
    })
    expected_df = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr1"],
        "exonStarts": [100, 200, 100, 200],
        "exonEnds": [150, 250, 150, 250],
        "strand": ["+", "-", "+", "-"],
        "exontype": ["skipped", "unique", "skipped", "unique"],
        "exon_position": ["first", "last", "first", "last"],
        "base_editor_name": ["MockBE1", "MockBE2", "MockBE1", "MockBE2"],
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
        "chrom": ["chr1"],
        "exonStarts": [100],
        "exonEnds": [150],
        "strand": ["+"],
        "exontype": ["skipped"],
        "exon_position": ["first"],
        "acceptor_sgrna_target_sequence": [["CCC+ATCG", "CCC+TAGC"]],
        "donor_sgrna_target_sequence": [["CCC+GCTA", "CCC+TACG"]]
    })
    input_df2 = pd.DataFrame({
        "chrom": ["chr1"],
        "exonStarts": [200],
        "exonEnds": [250],
        "strand": ["-"],
        "exontype": ["unique"],
        "exon_position": ["last"],
        "acceptor_sgrna_target_sequence": [["CCC+CGTA", "CCC+GCAT"]],
        "donor_sgrna_target_sequence": [["CCC+AGCT", "CCC+CTAG"]]
    })

    input_dict = {
        "MockBE1": input_df1,
        "MockBE2": input_df2
    }

    expected_df = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"],
        "exonStarts": [100, 100, 100, 100, 200, 200, 200, 200],
        "exonEnds": [150, 150, 150, 150, 250, 250, 250, 250],
        "strand": ["+", "+", "+", "+", "-", "-", "-", "-"],
        "exontype": ["skipped", "skipped", "skipped", "skipped", "unique", "unique", "unique", "unique"],
        "exon_position": ["first", "first", "first", "first", "last", "last", "last", "last"],
        "base_editor_name": ["MockBE1", "MockBE1", "MockBE1", "MockBE1", "MockBE2", "MockBE2", "MockBE2", "MockBE2"],
        "sgrna_target_sequence": ["CCC+ATCG", "CCC+TAGC", "CCC+GCTA", "CCC+TACG", "CCC+CGTA", "CCC+GCAT", "CCC+AGCT", "CCC+CTAG"],
        "site_type": ["acceptor", "acceptor", "donor", "donor", "acceptor", "acceptor", "donor", "donor"]
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
        "exontype": ["skipped", "skipped", "skipped", "skipped"],
        "exon_position": ["first", "first", "first", "first"],
        "base_editor_name": ["MockBE1", "MockBE1", "MockBE2", "MockBE2"],
        "sgrna_target_sequence": ["CCC+ATCG", "CCC+TAGC", "CCC+GCTA", "CCC+TACG"],
        "site_type": ["acceptor", "acceptor", "donor", "donor"]
    })

    mock_base_editors = [
    BaseEditor(
        base_editor_name="MockBE1",
        pam_sequence="NGG",
        editing_window_start_in_grna=10,
        editing_window_end_in_grna=15,
        base_editor_type="cbe"
    ),
    BaseEditor(
        base_editor_name="MockBE2",
        pam_sequence="NGG",
        editing_window_start_in_grna=15,
        editing_window_end_in_grna=20,
        base_editor_type="abe"
    )]
    expected_df = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr1"],
        "exonStarts": [100, 100, 100, 100],
        "exonEnds": [150, 150, 150, 150],
        "strand": ["+", "+", "+", "+"],
        "exontype": ["skipped", "skipped", "skipped", "skipped"],
        "exon_position": ["first", "first", "first", "first"],
        "base_editor_name": ["MockBE1", "MockBE1", "MockBE2", "MockBE2"],
        "base_editor_pam": ["NGG", "NGG", "NGG", "NGG"],
        "base_editor_editing_window_start": [10, 10, 15, 15],
        "base_editor_editing_window_end": [15, 15, 20, 20],
        "base_editor_type": ["cbe", "cbe", "abe", "abe"]
        "sgrna_target_sequence": ["CCC+ATCG", "CCC+TAGC", "CCC+GCTA", "CCC+TACG"],
        "site_type": ["acceptor", "acceptor", "donor", "donor"],
    })
