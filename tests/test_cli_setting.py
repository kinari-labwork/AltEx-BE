import pytest
import argparse
from altex_be.cli_setting import parse_base_editors,check_input_output_directories, get_base_editors_from_args, is_supported_assembly_name_in_crispr_direct, parse_gene_file
from altex_be.sgrna_designer import BaseEditor

def test_parse_base_editors_valid():
    args = argparse.Namespace(
        base_editor_name="TestBE",
        base_editor_pam="NGG",
        base_editor_window_start="10",
        base_editor_window_end="15",
        base_editor_type="cbe"
    )
    result = parse_base_editors(args)
    assert isinstance(result, dict)
    assert "TestBE" in result
    be = result["TestBE"]
    assert isinstance(be, BaseEditor)
    assert be.base_editor_name == "TestBE"
    assert be.pam_sequence == "NGG"
    assert be.editing_window_start_in_grna == 10
    assert be.editing_window_end_in_grna == 15
    assert be.base_editor_type == "cbe"

def test_parse_gene_file(tmp_path):
    # ダミーの遺伝子ファイル作成
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("GeneA\nGeneB\n  \nGeneC\n")
    result = parse_gene_file(gene_file)
    assert isinstance(result, set)
    assert "GeneA" in result
    assert "GeneB" in result
    assert "GeneC" in result
    assert len(result) == 3  # 空行が除外されていることを確認

def test_get_base_editors_from_args(tmp_path):
    # ダミーのbase editorファイル作成
    be_file = tmp_path / "editors.csv"
    be_file.write_text(
        "base_editor_name,pam_sequence,editing_window_start,editing_window_end,base_editor_type\n"
        "TestBE,NGG,10,15,cbe\n"
    )
    args = argparse.Namespace(be_files=str(be_file))
    result = get_base_editors_from_args(args)
    assert isinstance(result, dict)
    assert "TestBE" in result
    be = result["TestBE"]
    assert isinstance(be, BaseEditor)
    assert be.base_editor_name == "TestBE"

def test_get_base_editors_from_args_invalid(tmp_path):
    be_file = tmp_path / "invalid_editors.csv"
    # pam sequence 列が欠けているパターン
    be_file.write_text(
        "base_editor_name,editing_window_start,editing_window_end,base_editor_type\n"
        "TestBE,10,15,cbe\n"
    )
    args = argparse.Namespace(be_files=str(be_file))
    with pytest.raises(ValueError):
        get_base_editors_from_args(args)

def test_check_input_output_directories(tmp_path):
    # 1. Create dummy input files and an output directory
    refflat_file = tmp_path / "refflat.txt"
    refflat_file.touch()
    fasta_file = tmp_path / "genome.fa"
    fasta_file.touch()
    output_dir = tmp_path / "output"

    output_dir.mkdir()

    # 2. Test with valid directories
    check_input_output_directories(refflat_file, fasta_file, output_dir)
    assert output_dir.is_dir()

    # 3. Test with a non-existent input file
    with pytest.raises(FileNotFoundError):
        check_input_output_directories(
            tmp_path / "non_existent.txt", fasta_file, output_dir
        )

def test_is_supported_assembly_name_in_crispr_direct():
    assert is_supported_assembly_name_in_crispr_direct("TestAssembly") is False
    assert is_supported_assembly_name_in_crispr_direct("mm39") is True
