import pandas as pd
from altex_be.bed_for_ucsc_custom_track_maker import format_sgrna_for_ucsc_custom_track

def test_format_sgrna_for_ucsc_custom_track():
    # 1. Create a sample input DataFrame
    data = {
        'geneName': ['MYGENE'],
        'site_type': ['donor'],
        'base_editor_name': ['ABE8e'],
        'base_editor_type': ['abe'],
        'uuid': ['1234-5678'],
        'sgrna_sequence': ['ATGCATGCATGCATGCATGC'],
        'pam+20bp_exact_match_count': [1],
        'strand': ['+'],
        'sgrna_start_in_genome': [100],
        'sgrna_end_in_genome': [120],
        'chrom': ['chr1']
    }
    input_df = pd.DataFrame(data)

    # 2. Call the function with the sample data
    output_df = format_sgrna_for_ucsc_custom_track(input_df)

    # 3. Assert that the output is correct
    assert 'chrom' in output_df.columns
    assert 'chromStart' in output_df.columns
    assert 'chromEnd' in output_df.columns
    assert 'name' in output_df.columns
    assert 'score' in output_df.columns
    assert 'strand' in output_df.columns
    assert output_df['chromStart'].iloc[0] == 100
    assert output_df['chromEnd'].iloc[0] == 120
    assert output_df['name'].iloc[0] == 'MYGENE_donor_ABE8e_1234-5678'
    assert output_df['score'].iloc[0] == 1
    assert output_df['strand'].iloc[0] == '+'

def test_format_sgrna_for_ucsc_custom_track_edge_case():
    # 1. Create a sample input DataFrame
    data = {
        'geneName': ['MYGENE'],
        'site_type': ['donor'],
        'base_editor_name': ['ABE8e'],
        'base_editor_type': ['abe'],
        'uuid': ['1234-5678'],
        'sgrna_sequence': ['ATGCATGCATGCATGCATGC'],
        'pam+20bp_exact_match_count': [101],
        'strand': ['+'],
        'sgrna_start_in_genome': [100],
        'sgrna_end_in_genome': [120],
        'chrom': ['chr1']
    }
    input_df = pd.DataFrame(data)

    # 2. Call the function with the sample data
    output_df = format_sgrna_for_ucsc_custom_track(input_df)

    # 3. Assert that the output is correct
    assert 'chrom' in output_df.columns
    assert 'chromStart' in output_df.columns
    assert 'chromEnd' in output_df.columns
    assert 'name' in output_df.columns
    assert 'score' in output_df.columns
    assert 'strand' in output_df.columns
    assert output_df['chromStart'].iloc[0] == 100
    assert output_df['chromEnd'].iloc[0] == 120
    assert output_df['name'].iloc[0] == 'MYGENE_donor_ABE8e_1234-5678'
    assert output_df['score'].iloc[0] == 100
    assert output_df['strand'].iloc[0] == '+'