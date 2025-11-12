from dataclasses import (dataclass)

@dataclass(frozen=True)
class BaseEditor:
    """
    BaseEditorの情報を保持するためのdataclass
    """
    base_editor_name: str # BaseEditorの名前
    pam_sequence: str # PAM配列
    editing_window_start_in_grna: int # 編集ウィンドウの開始位置 (1-indexed)
    editing_window_end_in_grna: int # 編集ウィンドウの終了位置 (1-indexed)
    base_editor_type: str # "CBE" or "ABE"

PRESET_BASE_EDITORS = {
        "target_aid_ngg": BaseEditor(
            base_editor_name="target_aid_ngg",
            pam_sequence="NGG",
            editing_window_start_in_grna=17,
            editing_window_end_in_grna=19,
            base_editor_type="cbe"
        ),
        "target_aid_ng": BaseEditor(
            base_editor_name="target_aid_ng",
            pam_sequence="NG",
            editing_window_start_in_grna=17,
            editing_window_end_in_grna=19,
            base_editor_type="cbe"
        ),
        "be4max_ngg": BaseEditor(
            base_editor_name="be4max_ngg",
            pam_sequence="NGG",
            editing_window_start_in_grna=12,
            editing_window_end_in_grna=17,
            base_editor_type="cbe"
        ),
        "be4max_ng": BaseEditor(
            base_editor_name="be4max_ng",
            pam_sequence="NG",
            editing_window_start_in_grna=12,
            editing_window_end_in_grna=17,
            base_editor_type="cbe"
        ),
        "abe8e_ngg": BaseEditor(
            base_editor_name="abe8e_ngg",
            pam_sequence="NGG",
            editing_window_start_in_grna=12,
            editing_window_end_in_grna=17,
            base_editor_type="abe"
        ),
        "abe8e_ng": BaseEditor(
            base_editor_name="abe8e_ng",
            pam_sequence="NG",
            editing_window_start_in_grna=12,
            editing_window_end_in_grna=17,
            base_editor_type="abe"
        ),
    }