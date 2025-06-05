import pandas as pd
#export PYTHONPATH="$PYTHONPATH:./src" をbashで実行すること
from preprocess_and_additional_annotation import (
    drop_abnormal_mapped_transcripts,
    cording_information_annotator,
    flame_information_annotator,
    max_min_exon_count_annotator
)

#export PYTHONPATH="$PYTHONPATH:./src" をbashで実行すること

data = pd.read_csv("data/exon_classification.csv")

#transcript nameが重複している列を削除
data = data.drop_duplicates(subset=["name"],keep=False)
print(data.shape)

# 異常な染色体にマッピングされたトランスクリプトを削除
data = drop_abnormal_mapped_transcripts(data)

# コーディング情報を追加
data = cording_information_annotator(data)
# フレーム情報を追加
data = flame_information_annotator(data)
# 最大および最小のエクソン数を追加
data = max_min_exon_count_annotator(data)
# データフレームを保存
data.to_csv("data/exon_classification_with_additional_info.csv", index=False)