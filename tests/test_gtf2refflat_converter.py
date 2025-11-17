import pandas as pd
from pathlib import Path
from altex_be.gtf2refflat_converter import (
    gtf_to_refflat,
)

def test_gtf_to_refflat():
    assembly_name = "mm39"
    dummy_gtf_path = Path("tests/data/test.gtf")
    output_path = Path("tests/data")
    gtf_to_refflat(dummy_gtf_path, output_path, assembly_name=assembly_name)
    refflat_df = pd.read_csv(output_path / f"converted_refflat_{assembly_name}.txt", sep="\t", header=None)
    assert not refflat_df.empty
    assert refflat_df.shape[1] == 11  # geneName列が追加されていることを確認
    assert refflat_df.iloc[0,0].startswith("Gm")  # geneName列に正しい遺伝子記号が追加されていることを確認
    (output_path / f"converted_refflat_{assembly_name}.txt").unlink()  # クリーンアップ
    return