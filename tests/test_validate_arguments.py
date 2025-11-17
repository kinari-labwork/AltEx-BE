import pytest
import argparse
from altex_be.manage_arguments.validate_arguments import (
    is_supported_assembly_name_in_crispr_direct,
    is_input_output_directories,
    is_base_editors_provided,
    is_interest_genes_provided,
)
from altex_be.sgrna_designer import BaseEditor

def test_check_input_output_directories(tmp_path):
    # 1. Create dummy input files and an output directory
    refflat_file = tmp_path / "refflat.txt"
    refflat_file.touch()
    gtf_file = tmp_path / "annotations.gtf"
    gtf_file.touch()
    fasta_file = tmp_path / "genome.fa"
    fasta_file.touch()
    output_dir = tmp_path / "output"

    output_dir.mkdir()
    parser = argparse.ArgumentParser()
    # 2. Test with valid directories
    is_input_output_directories(refflat_file, None, fasta_file, output_dir, parser)
    is_input_output_directories(None, gtf_file, fasta_file, output_dir, parser)
    assert output_dir.is_dir()

    # 3. Test with a non-existent input file
    with pytest.raises(SystemExit):
        is_input_output_directories(
            tmp_path / "non_existent.txt", None, fasta_file, output_dir, parser
        )

def test_supported_assembly(caplog):
    # サポートされているアセンブリ名
    assert is_supported_assembly_name_in_crispr_direct("mm39") is None
    # ログに警告が出ていないことを確認
    assert "not supported by CRISPRdirect" not in caplog.text

def test_unsupported_assembly(caplog):
    # サポートされていないアセンブリ名
    with caplog.at_level("WARNING"):
        assert is_supported_assembly_name_in_crispr_direct("unknown_assembly") is None
        # 警告ログが出ていることを確認
        assert "unknown_assembly is not supported by CRISPRdirect" in caplog.text

def test_is_base_editors_provided():
    parser = argparse.ArgumentParser()
    # 1. Test with empty base editors
    with pytest.raises(SystemExit):
        is_base_editors_provided({}, parser)

    # 2. Test with non-empty base editors
    base_editors = {"BE1": BaseEditor(
        base_editor_name="BE1",
        pam_sequence="NGG",
        editing_window_start_in_grna=10,
        editing_window_end_in_grna=15,
        base_editor_type="cbe"
    )}
    # Should not raise an error
    is_base_editors_provided(base_editors, parser)

def test_is_interest_genes_provided():
    parser = argparse.ArgumentParser()
    # 1. Test with empty interest gene list
    with pytest.raises(SystemExit):
        is_interest_genes_provided([], parser)

    # 2. Test with non-empty interest gene list
    interest_genes = ["BRCA1", "TP53"]
    # Should not raise an error
    is_interest_genes_provided(interest_genes, parser)

