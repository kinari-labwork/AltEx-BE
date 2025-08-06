from altex_aid.sgrna_designer import BaseEditor

def parse_base_editors(base_editor_args:list[str], base_editors: list) -> list[BaseEditor]: 
    for editor_arg in base_editor_args:
        try:
            name, pam, window_start, window_end, editor_type = editor_arg.split(",")
            base_editors.append(
                BaseEditor(
                    base_editor_name=name,
                    pam_sequence=pam,
                    editing_window_start_in_grna=int(window_start),
                    editing_window_end_in_grna=int(window_end),
                    base_editor_type=editor_type.lower()
                )
            )
        except ValueError:
            print(
                f"Invalid BaseEditor format: {editor_arg}. Skipping."
                "Expected format: str[name],str[pam],int[window_start],int[window_end],str[editor_type]"
            )

    return base_editors