from altex_aid.sgrna_designer import BaseEditor
import argparse

def parse_base_editors(args: argparse.Namespace) -> list[BaseEditor] | None:
    base_editors = []
    if all([args.be_n, args.be_p, args.be_ws, args.be_we, args.be_t]):
        try:
            base_editors.append(
                BaseEditor(
                    base_editor_name=args.be_n,
                    pam_sequence=args.be_p.upper(),
                    editing_window_start_in_grna=int(args.be_ws),
                    editing_window_end_in_grna=int(args.be_we),
                    base_editor_type=args.be_t.lower(),
            )
        )
            return base_editors
        except ValueError as e:
            print(f"Error parsing base editor information: {e}")
            return None
    else:
        print("Base editor information is incomplete. Please provide all required parameters.")
        return None

def show_base_editors_info(base_editors: list[BaseEditor]):
    if base_editors is None:
        print("No base editors available to display.")
        return

    for base_editor in base_editors:
        print(f"  - {base_editor.base_editor_name} (Type: {base_editor.base_editor_type}, PAM: {base_editor.pam_sequence}, "
            f"Window: {base_editor.editing_window_start_in_grna}-{base_editor.editing_window_end_in_grna})")