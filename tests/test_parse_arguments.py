import pytest
import argparse
from altex_be.manage_arguments.parse_arguments import (
    parse_path_from_args,
    parse_gene_file,
    parse_assembly_name_from_args,
    parse_genes_from_args,
    parse_base_editors_from_args,
    parse_base_editors_from_file,
    parse_base_editors_from_presets,
)
from altex_be.sgrna_designer import BaseEditor

def test_parse_path_from_args(tmp_path):
    refflat_file = tmp_path / "refflat.txt"
    fasta_file = tmp_path / "genome.fa"
    output_dir = tmp_path / "output"
    gtf_file = tmp_path / "annotations.gtf"

    args = argparse.Namespace(
        refflat_path=str(refflat_file),
        fasta_path=str(fasta_file),
        gtf_path=str(gtf_file),
        output_dir=str(output_dir)
    )

    refflat_path, gtf_path, fasta_path, output_directory = parse_path_from_args(args)
    assert refflat_path == refflat_file
    assert gtf_path == gtf_file
    assert fasta_path == fasta_file
    assert output_directory == output_dir

def test_parse_assembly_name_from_args():
    args = argparse.Namespace(assembly_name="hg38")
    assembly_name = parse_assembly_name_from_args(args)
    assert assembly_name == "hg38"


def test_parse_base_editors_from_args_valid():
    args = argparse.Namespace(
        be_name="TestBE",
        be_pam="NGG",
        be_window_start="10",
        be_window_end="15",
        be_type="cbe"
    )
    base_editors = {}
    parser = argparse.ArgumentParser()
    result = parse_base_editors_from_args(args, parser, base_editors)
    assert isinstance(result, dict)
    assert "TestBE" in result
    be = result["TestBE"]
    assert isinstance(be, BaseEditor)
    assert be.base_editor_name == "TestBE"
    assert be.pam_sequence == "NGG"
    assert be.editing_window_start_in_grna == 10
    assert be.editing_window_end_in_grna == 15
    assert be.base_editor_type == "cbe"

def test_parse_base_editors_from_args_invalid():
    args = argparse.Namespace(
        be_name="TestBE",
        be_pam=None,
        be_window_start="10",
        be_window_end="15",
        be_type="cbe"
    )
    base_editors = {}
    parser = argparse.ArgumentParser()
    with pytest.raises(SystemExit):
        parse_base_editors_from_args(args, parser, base_editors)

def test_parse_base_editors_from_presets():
    args = argparse.Namespace(be_preset="abe8e_ngg")
    parser = argparse.ArgumentParser()
    result = parse_base_editors_from_presets(args, parser, {})
    assert isinstance(result, dict)
    print(result)
    be = result["abe8e_ngg"]
    assert isinstance(be, BaseEditor)
    assert be.base_editor_name == "abe8e_ngg"

def test_parse_gene_file(tmp_path):
    # ダミーの遺伝子ファイル作成
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("GeneA\nGeneB\n  \nGeneC\n")
    result = parse_gene_file(gene_file)
    assert isinstance(result, list)
    assert "GeneA" in result
    assert "GeneB" in result
    assert "GeneC" in result
    assert len(result) == 3  # 空行が除外されていることを確認

def test_parse_genes_from_args():
    # ダミーの遺伝子ファイル作成
    args = argparse.Namespace(
        gene_symbols=["GeneA", "GeneB"],
        refseq_ids=["NM_001", "NM_002"],
        ensembl_ids=["ENST000003", "ENST000004"],
        gene_file=None
    )
    parser = argparse.ArgumentParser()
    result = parse_genes_from_args(args, parser)
    assert isinstance(result, list)
    assert "GeneA" in result
    assert "GeneB" in result
    assert "NM_001" in result
    assert "NM_002" in result
    assert "ENST000003" in result
    assert "ENST000004" in result
    assert len(result) == 6

def test_parse_base_editors_from_files(tmp_path):
    # ダミーのbase editorファイル作成
    be_file = tmp_path / "editors.csv"
    be_file.write_text(
        "base_editor_name,pam_sequence,editing_window_start,editing_window_end,base_editor_type\n"
        "TestBE,NGG,10,15,cbe\n"
    )
    parser = argparse.ArgumentParser()
    args = argparse.Namespace(base_editor_files=str(be_file))
    result = parse_base_editors_from_file(args, parser)
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
    args = argparse.Namespace(base_editor_files=str(be_file))
    parser = argparse.ArgumentParser()
    with pytest.raises(SystemExit):
        parse_base_editors_from_file(args, parser)
