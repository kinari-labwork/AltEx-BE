import pandas as pd
#export PYTHONPATH="$PYTHONPATH:./src" をbashで実行すること
from altex_aid.refflat_preprocessor import (
    parse_exon_coordinates,
    calculate_exon_lengths,
    drop_abnormal_mapped_transcripts,
    cording_information_annotator,
    flame_information_annotator,
    variant_count_annotator,
    add_exon_position_flags
)

annotation_genome = "mm39"  # ここは必要に応じて変更してください

data = pd.read_csv(
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
        "exonEnds"])


#transcript nameが重複している列を削除
data = data.drop_duplicates(subset=["name"],keep=False)

# "exonStarts" と　"exonEnds"を扱いやすいように (start, end) のリストに変換する
data = parse_exon_coordinates(data)

# 各エキソンの長さを計算して追加
data = calculate_exon_lengths(data)

# 異常な染色体にマッピングされたトランスクリプトを削除
data = drop_abnormal_mapped_transcripts(data)

# コーディング情報を追加
data = cording_information_annotator(data)

# フレーム情報を追加
data = flame_information_annotator(data)

#遺伝子ごとのバリアントの数を追加
data = variant_count_annotator(data)

# エキソンの位置を示す列を追加
data = add_exon_position_flags(data)

# データフレームを保存
data.to_pickle("data/processed_refflat.pkl")