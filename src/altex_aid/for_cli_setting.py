from altex_aid.sgrna_designer import BaseEditor

def input_base_editors(base_editors: list[BaseEditor]) -> list[BaseEditor]:
    print("you can input multiple BaseEditor information. To finish, just press Enter without typing anything.")
    while True:
        name = input("name of BaseEditor (or press Enter to finish): ")
        if not name:
            break
        pam = input("PAM sequence (e.g., NGG): ")
        window_start = input("Editing window start position, count from PAM (e.g., 4): ")
        window_end = input("Editing window end position count from PAM (e.g., 8): ")
        editor_type = input("BaseEditor type (cbe or abe): ")
        # 必須項目が空ならスキップ
        if not (pam and window_start and window_end and editor_type):
            print("All fields are required. Please try again.")
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