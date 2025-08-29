from altex_be.sgrna_designer import BaseEditor
import pandas as pd
import uuid
import logging
from . import logging_config # noqa: F401


def convert_empty_list_into_na(target_exon_with_sgrna_dict: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """
    Convert empty lists in the DataFrame to NaN.
    """
    for df in target_exon_with_sgrna_dict.values():
        for col in df.columns:
            df[col] = df[col].apply(lambda x: pd.NA if isinstance(x, list) and not x else x)
    return target_exon_with_sgrna_dict

def prepare_melted_df(target_exon_with_sgrna_df: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose : 1行-複数のsgRNAの状態からexplodeして1行-1sgRNAに変換する前処理として、列名を整形し、データフレームを整形する
    Parameters : target_exon_with_sgrna_df: あるBEに対して設計されたsgRNAのdf。 acceptorとdonorの2つのサイトタイプの情報が含まれる
    Return : melted_df: 整形後のデータフレーム
    """

    foundation_cols = [
        "geneName",
        "chrom",
        "exonStarts",
        "exonEnds",
        "strand",
        "exonlengths",
        "coding",
        "flame",
        "exontype",
        "exon_position",
        "uuid"
        ]
    df_list = []
    for site in ["acceptor", "donor"]:
        cols = [col for col in target_exon_with_sgrna_df.columns if col.startswith(site)]
        site_df = target_exon_with_sgrna_df[foundation_cols + cols].rename(
            columns={
                col: f"{col.replace(f'{site}_', '')}" for col in cols
            }
        )
        site_df = site_df.dropna(subset=["sgrna_target_sequence"])
        site_df["site_type"] = site
        df_list.append(site_df.reset_index(drop=True))
    melted_df = pd.concat(df_list, ignore_index=True)
    return melted_df

def is_sgrna_designed(melted_df: pd.DataFrame) -> bool:
    """
    Purpose: Validate the exploded sgRNA DataFrame.
    """
    if melted_df.empty:
        logging.warning("there are no designed sgRNAs for your interest gene and your base editor.")
        return False
    return True

def explode_sgrna_df(target_exon_with_sgrna_dict: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Purpose: 複数のBEに対して設計されたsgRNAのdfが、1列-1sgRNAになるように変換する
    """
    # 各BEに対して、sgRNAのdfを1列-1sgRNAに変換
    exploded_dfs = []
    for be, df in target_exon_with_sgrna_dict.items():
        melted_df = prepare_melted_df(df)
        melted_df["base_editor_name"] = be
        sgrna_cols = [col for col in melted_df.columns if col.startswith("sgrna_")]
        melted_df = melted_df.explode(column=sgrna_cols)
        exploded_dfs.append(melted_df)
    # すべてのBEのdfを結合
    exploded_sgrna_df = pd.concat(exploded_dfs, ignore_index=True)
    return exploded_sgrna_df

def add_base_editor_info_to_df(exploded_sgrna_df: pd.DataFrame, base_editors: dict[str, BaseEditor]) -> pd.DataFrame:
    """
    Purpose: exploded_sgrna_dfにBaseEditorの情報を追加する
    """
    # BaseEditorのリストをDataFrameに変換
    be_info = [
        {
            "base_editor_name": be.base_editor_name,
            "base_editor_pam_sequence": be.pam_sequence,
            "base_editor_editing_window_start": be.editing_window_start_in_grna,
            "base_editor_editing_window_end": be.editing_window_end_in_grna,
            "base_editor_type": be.base_editor_type,
        }
        for be in base_editors
    ]
    be_df = pd.DataFrame(be_info)

    # 'base_editor_name'をキーとしてマージ（結合）
    exploded_sgrna_df = pd.merge(exploded_sgrna_df, be_df, on="base_editor_name", how="left")
    return exploded_sgrna_df

def update_uuid_unique_to_every_sgrna(exploded_sgrna_df: pd.DataFrame) -> pd.DataFrame:
    """
    Purpose: ここまでの処理は、情報のマージのために各エキソンに対してユニークなUUIDを設定していたが、ここからはsgRNAごとにユニークなUUIDを設定する
    """
    exploded_sgrna_df['uuid'] = [uuid.uuid4().hex for _ in range(len(exploded_sgrna_df))]
    return exploded_sgrna_df

def get_pam_sequence_based_on_real_sequence(exploded_sgrna_df:pd.DataFrame) -> pd.DataFrame:
    """
    Purpose: 今の時点では、base_editor_pam_sequence列は、"NGG"のように曖昧配列で表現されている。
            これを、実際のsgRNA配列に基づいて、"AGG"や"TGG"のように具体的な配列の表現の列を追加する
    """
    exploded_sgrna_df = exploded_sgrna_df.reset_index(drop=True)
    pam_seq_list = []

    for index, row in exploded_sgrna_df.iterrows():
        if row["site_type"] == "acceptor" and row["base_editor_type"] == "abe":
            pam_seq = row["sgrna_sequence"].split("+")[0]
        else:
            pam_seq = row["sgrna_sequence"].split("+")[1]
        pam_seq_list.append(pam_seq)
    col_index = exploded_sgrna_df.columns.get_loc("sgrna_sequence") + 1
    exploded_sgrna_df.insert(col_index, 'pam_sequence', pam_seq_list)

def format_output(target_exon_with_sgrna_dict: dict[str, pd.DataFrame], 
                base_editors: dict[str, BaseEditor]) -> pd.DataFrame:
    """
    Purpose: このモジュールのラップ関数
    """
    # 空のリストをNaNに変換
    target_exon_with_sgrna_dict = convert_empty_list_into_na(target_exon_with_sgrna_dict)

    # sgRNAのdfをexplodeして1列-1sgRNAに変換
    exploded_sgrna_df = explode_sgrna_df(target_exon_with_sgrna_dict)

    # exploded_sgrna_dfの検証
    if not is_sgrna_designed(exploded_sgrna_df):
        return pd.DataFrame()  # 空のDataFrameを返す

    # BaseEditorの情報を追加
    exploded_sgrna_df = add_base_editor_info_to_df(exploded_sgrna_df, base_editors)

    # UUIDをsgRNAごとにユニークに更新
    exploded_sgrna_df = update_uuid_unique_to_every_sgrna(exploded_sgrna_df)

    # geneNameを基にソート
    exploded_sgrna_df = exploded_sgrna_df.sort_values(by=["geneName", "exon_position"])

    return exploded_sgrna_df.reset_index(drop=True)
