import pandas as pd
from pathlib import Path
from altex_be.gtf2refflat_converter import (
    run_gtf2genepred,
    convert_gtf_to_refflat_format,
    add_genesymbol_to_imcomplete_refflat
)

def test_run_gtf2genepred():
    dummy_gtf_path = Path("tests/data/test.gtf")
    dummy_genepred_path = Path("tests/data/test.genepred")
    run_gtf2genepred(dummy_gtf_path, dummy_genepred_path)
    assert dummy_genepred_path.with_suffix('.genepred').exists()
    genepred = pd.read_csv(dummy_genepred_path.with_suffix('.genepred'), sep="\t", header=None)
    assert not genepred.empty
    return

def test_convert_gtf_to_refflat_format():
    dummy_genepred_path = Path("tests/data/test.genepred")
    refflat_df_without_genesymbol = convert_gtf_to_refflat_format(dummy_genepred_path)
    assert not refflat_df_without_genesymbol.empty
    assert refflat_df_without_genesymbol.shape[1] == 10
    assert refflat_df_without_genesymbol.iloc[0, 1].startswith("chr")  # 染色体列に"chr"接頭辞が追加されていることを確認
    return

def test_add_genesymbol_to_imcomplete_refflat():
    refflat_without_genesymbol = convert_gtf_to_refflat_format(Path("tests/data/test.genepred"))
    dummy_gtf_path = Path("tests/data/test.gtf")
    refflat_with_genesymbol = add_genesymbol_to_imcomplete_refflat(refflat_without_genesymbol, dummy_gtf_path)
    assert not refflat_with_genesymbol.empty
    assert refflat_with_genesymbol.shape[1] == 11  # geneName列が追加されていることを確認
    assert refflat_with_genesymbol.iloc[0,0].startswith("Gm")  # geneName列に正しい遺伝子記号が追加されていることを確認
    return