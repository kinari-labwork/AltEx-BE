import pandas as pd

# export PYTHONPATH="$PYTHONPATH:./src" をbashで実行すること
from altex_aid.refflat_preprocessor import (
    add_exon_position_flags,
    annotate_cording_information,
    annotate_flame_information,
    calculate_exon_lengths,
    drop_abnormal_mapped_transcripts,
    parse_exon_coordinates,
)

annotation_genome = "mm39"  # ここは必要に応じて変更してください

refflat = pd.read_csv(
    f"data/{annotation_genome}/refFlat.txt",
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


# transcript nameが重複している列を削除
refflat = refflat.drop_duplicates(subset=["name"], keep=False)

# "exonStarts" と　"exonEnds"を扱いやすいように (start, end) のリストに変換する
refflat = parse_exon_coordinates(refflat)

# 各エキソンの長さを計算して追加
refflat = calculate_exon_lengths(refflat)

# 異常な染色体にマッピングされたトランスクリプトを削除
refflat = drop_abnormal_mapped_transcripts(refflat)

# コーディング情報を追加
refflat = annotate_cording_information(refflat)

# フレーム情報を追加
refflat = annotate_flame_information(refflat)

# エキソンの位置を示す列を追加
refflat = add_exon_position_flags(refflat)

# データフレームを保存
refflat.to_pickle("data/processed_refflat.pkl")
