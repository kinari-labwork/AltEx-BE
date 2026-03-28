import pandas as pd
import numpy as np
import pytest
from altex_be.sgrna_prioritizer import prioritize_sgrna


@pytest.fixture
def sample_sgrna_df():
    """Create a sample sgRNA DataFrame for testing."""
    return pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene1", "gene2", "gene2"],
        "exonStarts": [100, 100, 100, 200, 200],
        "exonEnds": [150, 150, 150, 250, 250],
        "sgrna_sequence": ["ATGCATGCATGCATGCATGC", "ATGCATGCATGCATGCATGC", "GGGGGGGGGGGGGGGGGGGG", "ATGCATGCATGC", "CCCCCCCCCCCCCCCCCCCC"],
        "pam+20bp_exact_match_count": [1, 5, 1, 2, 1],
        "pam+12bp_exact_match_count": [10, 20, 10, 15, 10],
        "sgrna_overlap_between_cds_and_editing_window": [0, 0, 1, 0, 0],
        "sgrna_possible_unintended_edited_base_count": [0, 0, 2, 0, 0],
    })


def test_prioritize_sgrna_basic(sample_sgrna_df):
    """Test basic prioritization with sample data."""
    result = prioritize_sgrna(sample_sgrna_df)
    
    # Check that output has required columns
    assert "sgrna_priority" in result.columns
    
    # Check that intermediate columns were removed
    assert "gc_valid" not in result.columns
    assert "gc_percentage" not in result.columns
    assert "composite_score" not in result.columns
    
    # Check shape is preserved
    assert len(result) == len(sample_sgrna_df)


def test_prioritization_order_by_20bp_offtarget(sample_sgrna_df):
    """Test that sgRNAs are primarily prioritized by 20bp off-target count."""
    result = prioritize_sgrna(sample_sgrna_df)
    
    # Within each exon group (gene1), sgRNA with off-target=1 should rank higher than off-target=5
    gene1_results = result[result["geneName"] == "gene1"].sort_values("sgrna_priority")
    
    # First sgRNA (off-target 20bp = 1, GC valid, no CDS) should be rank 1
    assert gene1_results.iloc[0]["pam+20bp_exact_match_count"] == 1
    assert gene1_results.iloc[0]["sgrna_priority"] == 1


def test_gc_content_validation():
    """Test that GC content is correctly validated (40-60% range)."""
    # GC% of sequence: ATGCATGCATGCATGCATGC (10G + 10C out of 20 = 50%)
    gc_valid_seq = "ATGCATGCATGCATGCATGC"
    
    # GC% of sequence: AAAAAAAAAAAAAAAAAAAA (0G + 0C out of 20 = 0%)
    gc_invalid_seq_low = "AAAAAAAAAAAAAAAAAAAA"
    
    # GC% of sequence: GGGGGGGGGGGGGGGGGGGG (20G + 0C out of 20 = 100%)
    gc_invalid_seq_high = "GGGGGGGGGGGGGGGGGGGG"
    
    df = pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene1"],
        "exonStarts": [100, 100, 100],
        "exonEnds": [150, 150, 150],
        "sgrna_sequence": [gc_valid_seq, gc_invalid_seq_low, gc_invalid_seq_high],
        "pam+20bp_exact_match_count": [1, 1, 1],  # All same to test secondary criterion
        "pam+12bp_exact_match_count": [10, 10, 10],
        "sgrna_overlap_between_cds_and_editing_window": [0, 0, 0],
        "sgrna_possible_unintended_edited_base_count": [0, 0, 0],
    })
    
    result = prioritize_sgrna(df)
    
    # GC valid (50%) should be ranked first
    assert result.iloc[0]["sgrna_sequence"] == gc_valid_seq
    assert result.iloc[0]["sgrna_priority"] == 1


def test_tiebreaker_12bp_offtarget():
    """Test that 12bp off-target is used as tiebreaker when 20bp is equal."""
    df = pd.DataFrame({
        "geneName": ["gene1", "gene1"],
        "exonStarts": [100, 100],
        "exonEnds": [150, 150],
        "sgrna_sequence": ["ATGCATGCATGCATGCATGC", "ATGCATGCATGCATGCATGC"],
        "pam+20bp_exact_match_count": [1, 1],  # Same
        "pam+12bp_exact_match_count": [10, 5],  # Different - 5 should rank higher
        "sgrna_overlap_between_cds_and_editing_window": [0, 0],
        "sgrna_possible_unintended_edited_base_count": [0, 0],
    })
    
    result = prioritize_sgrna(df)
    
    # sgRNA with 12bp count of 5 should be ranked first
    top_ranked = result.iloc[0]
    assert top_ranked["pam+12bp_exact_match_count"] == 5
    assert top_ranked["sgrna_priority"] == 1


def test_final_tiebreaker_unintended_edits():
    """Test that unintended edits are used as final tiebreaker."""
    df = pd.DataFrame({
        "geneName": ["gene1", "gene1"],
        "exonStarts": [100, 100],
        "exonEnds": [150, 150],
        "sgrna_sequence": ["ATGCATGCATGCATGCATGC", "ATGCATGCATGCATGCATGC"],
        "pam+20bp_exact_match_count": [1, 1],  # Same
        "pam+12bp_exact_match_count": [10, 10],  # Same
        "sgrna_overlap_between_cds_and_editing_window": [0, 0],
        "sgrna_possible_unintended_edited_base_count": [2, 0],  # Different - 0 should rank higher
    })
    
    result = prioritize_sgrna(df)
    
    # sgRNA with unintended edits = 0 should be ranked first
    top_ranked = result.iloc[0]
    assert top_ranked["sgrna_possible_unintended_edited_base_count"] == 0
    assert top_ranked["sgrna_priority"] == 1


def test_priority_reset_per_exon():
    """Test that priority ranking resets for each exon."""
    df = pd.DataFrame({
        "geneName": ["gene1", "gene1", "gene1", "gene2", "gene2"],
        "exonStarts": [100, 100, 100, 200, 200],
        "exonEnds": [150, 150, 150, 250, 250],
        "sgrna_sequence": ["ATGCATGCATGCATGCATGC"] * 5,
        "pam+20bp_exact_match_count": [3, 1, 2, 1, 2],
        "pam+12bp_exact_match_count": [10] * 5,
        "sgrna_overlap_between_cds_and_editing_window": [0] * 5,
        "sgrna_possible_unintended_edited_base_count": [0] * 5,
    })
    
    result = prioritize_sgrna(df)
    
    # Each exon should have its own priority ranking
    gene1_priorities = result[result["geneName"] == "gene1"]["sgrna_priority"].tolist()
    gene2_priorities = result[result["geneName"] == "gene2"]["sgrna_priority"].tolist()
    
    assert sorted(gene1_priorities) == [1, 2, 3]
    assert sorted(gene2_priorities) == [1, 2]


def test_composite_score_internally_calculated(sample_sgrna_df):
    """Test that composite score is used for sorting (internal calculation)."""
    result = prioritize_sgrna(sample_sgrna_df)
    
    # Test that sgRNA is prioritized correctly based on internal scoring
    gene1_results = result[result["geneName"] == "gene1"]
    
    # First result should have lowest 20bp count (1 instead of 5)
    assert gene1_results.iloc[0]["pam+20bp_exact_match_count"] == 1
