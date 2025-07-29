from altex_aid.sgrna_designer import BaseEditor

def input_base_editors():
    base_editors = []
    print("BaseEditor情報を複数入力できます。終了したい場合は何も入力せずEnterしてください。")
    while True:
        name = input("BaseEditor名（例: CBE, ABE）: ")
        if not name:
            break
        pam = input("PAM配列（例: NGG）: ")
        window_start = input("編集ウィンドウ開始位置（例: 4）: ")
        window_end = input("編集ウィンドウ終了位置（例: 8）: ")
        editor_type = input("BaseEditorタイプ（CBE or ABE）: ")
        # 必須項目が空ならスキップ
        if not (pam and window_start and window_end and editor_type):
            print("必要な情報が不足しています。もう一度入力してください。")
            continue
        base_editors.append(
            BaseEditor(
                base_editor_name=name,
                pam_sequence=pam,
                editing_window_start_in_grna=int(window_start),
                editing_window_end_in_grna=int(window_end),
                base_editor_type=editor_type.upper()
            )
        )
    return base_editors