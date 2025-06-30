import pandas as pd
import pickle
#export PYTHONPATH="$PYTHONPATH:./src" をbashで実行すること
from altex_aid.preprocess_and_additional_annotation import (
    drop_abnormal_mapped_transcripts,
    cording_information_annotator,
    flame_information_annotator,
    variant_count_annotator,
    add_exon_position_flags
)

# pklで読み込むと、evalを使わずにリストを読み込むことができる
data = pd.read_pickle("data/exon_classification.pkl")

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
data.to_pickle("data/exon_classification_with_additional_info.pkl")