import os
import shutil
import subprocess
import tempfile
import time
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
from typing import Optional, List

import pandas as pd
import streamlit as st

# 動的にインポートを試みる
try:
    from altex_be.class_def.base_editors import PRESET_BASE_EDITORS
except ImportError:
    # ローカルでの開発中など、パスが通っていない場合を考慮
    from class_def.base_editors import PRESET_BASE_EDITORS


# ========= App config =========
APP_TITLE = "AltEx-BE UI"
st.set_page_config(page_title=APP_TITLE, page_icon="🧬", layout="wide")

# ========= Helpers (tk dialogs) =========
def select_folder():
    """Open a dialog to select a folder."""
    root = tk.Tk()
    root.withdraw()
    folder_path = filedialog.askdirectory(master=root, title="Select output directory")
    root.destroy()
    return folder_path

def select_file(file_type: str = "fasta"):
    """Open a dialog to select a file."""
    root = tk.Tk()
    root.withdraw()
    if file_type == "fasta":
        filetypes = [("FASTA files", "*.fa *.fasta *.fna"), ("All files", "*.* ")]
    elif file_type == "annotation":
        filetypes = [("Annotation files", "*.gtf *.txt"), ("All files", "*.* ")]
    else:
        filetypes = [("All files", "*.* ")]
    
    file_path = filedialog.askopenfilename(
        master=root,
        filetypes=filetypes,
        title=f"Select {file_type.capitalize()} file",
    )
    root.destroy()
    return file_path

def is_writable_dir(path: str) -> bool:
    """Check if a path is a writable directory."""
    return os.path.isdir(path) and os.access(path, os.W_OK)

def is_valid_fasta(path: str) -> bool:
    """Check if a file has a valid FASTA extension."""
    p = path.lower()
    return p.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"))

def is_valid_annotation(path: str) -> bool:
    """Check if a file has a valid annotation extension."""
    p = path.lower()
    return p.endswith((".gtf", ".txt"))

def build_command(
    fasta_path: str,
    annotation_path: str,
    outdir: str,
    assembly_name: str,
    gene_symbols: Optional[List[str]] = None,
    gene_file: Optional[str] = None,
    be_name: Optional[str] = None,
    be_type: Optional[str] = None,
    be_pam: Optional[str] = None,
    be_start: Optional[int] = None,
    be_end: Optional[int] = None,
) -> list[str]:
    """Build the altex-be command from UI inputs."""
    cmd = ["altex-be"]

    # Input files
    if annotation_path.lower().endswith(".gtf"):
        cmd.extend(["--gtf-path", annotation_path])
    else:
        cmd.extend(["--refflat-path", annotation_path])
    
    cmd.extend(["--fasta-path", fasta_path])
    cmd.extend(["--output-dir", outdir])
    cmd.extend(["--assembly-name", assembly_name])

    # Gene inputs
    if gene_symbols:
        cmd.extend(["--gene-symbols", *gene_symbols])
    elif gene_file:
        cmd.extend(["--gene-file", gene_file])

    # Custom base editor
    if be_name and be_type and be_pam and be_start and be_end:
        cmd.extend([
            "--be-name", be_name,
            "--be-type", be_type,
            "--be-pam", be_pam,
            "--be-start", str(be_start),
            "--be-end", str(be_end),
        ])
        
    return cmd

# ========= Session state defaults =========
if "input_fasta" not in st.session_state:
    st.session_state.input_fasta = ""
if "input_annotation" not in st.session_state:
    st.session_state.input_annotation = ""
if "output_dir" not in st.session_state:
    st.session_state.output_dir = ""
if "last_cmd" not in st.session_state:
    st.session_state.last_cmd = []
if "log_text" not in st.session_state:
    st.session_state.log_text = ""
if "last_run_outdir" not in st.session_state:
    st.session_state.last_run_outdir = ""


# ========= Header =========
st.markdown(f"# {APP_TITLE}")
st.caption("Local UI wrapper for AltEx-BE CLI — configure arguments, run jobs, preview results.")

with st.container(border=True):
    st.markdown("### 📥 1. Input required files")
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("**FASTA**")
        if st.button("Select FASTA path", use_container_width=True):
            fp = select_file("fasta")
            if fp:
                st.session_state.input_fasta = fp
        st.session_state.input_fasta = st.text_input("FASTA path", value=st.session_state.input_fasta, help="Path to the input FASTA file.")

    with col2:
        st.markdown("**Transcript annotation**")
        if st.button("Select annotation path", use_container_width=True):
            ap = select_file("annotation")
            if ap:
                st.session_state.input_annotation = ap
        st.session_state.input_annotation = st.text_input("Annotation path", value=st.session_state.input_annotation, help="Path to GTF or RefFlat file.")

    with col3:
        st.markdown("**Output directory**")
        if st.button("Select output directory", use_container_width=True):
            od = select_folder()
            if od:
                st.session_state.output_dir = od
        st.session_state.output_dir = st.text_input("Output directory", value=st.session_state.output_dir, help="Path to the output directory.")

with st.container(border=True):
    st.markdown("### 🧬 2. Select Target Gene and assembly")
    target_species = st.text_input("Enter the assembly name (e.g., hg38, mm39)", help="Specify the assembly name.")
    gene_mode = st.radio("Select gene input method", ["Text input", "File upload"], index=0, horizontal=True)
    
    target_genes_text = None
    uploaded_gene_file = None

    if gene_mode == "Text input":
        target_genes_text = st.text_input("Enter target gene symbol(s) or ID(s) (comma-separated)", help="Specify one or more gene names separated by commas.")
    else: # File upload
        uploaded_gene_file = st.file_uploader("Upload a text file with gene symbols or IDs (one per line)", type=["txt", "csv", "tsv"])

with st.container(border=True):
    st.markdown("### ⚙️ 3. Select Base Editing tools")
    mode = st.radio("Select mode", ["Preset (default)", "Preset + Custom"], index=0, help="Use default editors or add a custom one.", horizontal=True)
    
    editor_name, pam_sequences, editing_type, window_start, window_end = None, None, None, None, None
    if mode == "Preset + Custom":
        c1, c2 = st.columns(2)
        with c1:
            editor_name = st.text_input("Custom Editor Name")
            pam_sequences = st.text_input("Custom PAM Sequence(s)")
        with c2:
            editing_type = st.selectbox("Editing Type", ["cbe", "abe"], index=0)
            window_start, window_end = st.slider("Editing Window", 1, 30, (13, 17), 1)

    with st.expander("🔍 Show preset base editors"):
        df = pd.DataFrame([vars(be) for be in PRESET_BASE_EDITORS.values()])
        st.dataframe(df, use_container_width=True, hide_index=True)

with st.container(border=True):
    st.markdown("### 🧪 4. Run AltEx-BE")
    run_clicked = st.button("▶ Run AltEx-BE", type="primary", use_container_width=True)
    st.caption("Logs will stream to the main panel below.")
    log_box = st.empty()

    def stream_process(cmd: list[str], cwd: str | None = None) -> int:
        st.session_state.log_text = f"$ {' '.join(cmd)}\n\n"
        st.session_state.run_state = "Running"
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, universal_newlines=True, cwd=cwd
        )
        for line in proc.stdout:
            st.session_state.log_text += line
            tail = st.session_state.log_text[-4000:]
            log_box.code(tail, language="bash")
        return proc.wait()

    # --- Input validation ---
    fasta_ok = bool(st.session_state.input_fasta) and os.path.exists(st.session_state.input_fasta) and is_valid_fasta(st.session_state.input_fasta)
    annotation_ok = bool(st.session_state.input_annotation) and os.path.exists(st.session_state.input_annotation) and is_valid_annotation(st.session_state.input_annotation)
    out_ok = bool(st.session_state.output_dir) and os.path.exists(st.session_state.output_dir) and is_writable_dir(st.session_state.output_dir)
    assembly_ok = bool(target_species)
    gene_input_ok = (gene_mode == "Text input" and bool(target_genes_text)) or \
                    (gene_mode == "File upload" and uploaded_gene_file is not None)

    if run_clicked:
        # --- Final validation before run ---
        if not fasta_ok: st.error("Please select a valid FASTA file.")
        elif not annotation_ok: st.error("Please select a valid annotation file (GTF/RefFlat).")
        elif not out_ok: st.error("Please select a writable output directory.")
        elif not assembly_ok: st.error("Please enter an assembly name.")
        elif not gene_input_ok: st.error("Please provide at least one target gene or a gene file.")
        else:
            ts = time.strftime("%Y%m%d_%H%M%S")
            run_outdir = str(Path(st.session_state.output_dir) / f"altexbe_{ts}")
            os.makedirs(run_outdir, exist_ok=True)
            st.session_state.last_run_outdir = run_outdir

            gene_list = None
            gene_file_path = None
            if gene_mode == "Text input":
                gene_list = [g.strip() for g in target_genes_text.split(",") if g.strip()]
            else: # File upload
                with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".txt", dir=run_outdir) as tmp:
                    content = uploaded_gene_file.read().decode("utf-8")
                    tmp.write(content)
                    gene_file_path = tmp.name

            cmd = build_command(
                fasta_path=st.session_state.input_fasta,
                annotation_path=st.session_state.input_annotation,
                outdir=run_outdir,
                assembly_name=target_species,
                gene_symbols=gene_list,
                gene_file=gene_file_path,
                be_name=editor_name if mode == "Preset + Custom" else None,
                be_type=editing_type if mode == "Preset + Custom" else None,
                be_pam=pam_sequences if mode == "Preset + Custom" else None,
                be_start=window_start if mode == "Preset + Custom" else None,
                be_end=window_end if mode == "Preset + Custom" else None,
            )
            st.session_state.last_cmd = cmd
            
            rc = stream_process(cmd)
            if rc == 0:
                st.toast("AltEx-BE completed successfully!", icon="✅")
            else:
                st.toast(f"AltEx-BE failed with exit code {rc}.", icon="❌")

# --- Tabs for status and results ---
tabs = st.tabs(["✅ Inputs Status", "📊 Results Preview", "📝 Full Log", "ℹ️ About"])

with tabs[0]:
    st.metric("FASTA", "OK" if fasta_ok else "Missing/Invalid")
    st.metric("Annotation", "OK" if annotation_ok else "Missing/Invalid")
    st.metric("Output Dir", "OK" if out_ok else "Missing/Not writable")
    st.metric("Assembly", "OK" if assembly_ok else "Missing")
    st.metric("Gene Input", "OK" if gene_input_ok else "Missing")

with tabs[1]:
    st.caption(f"Last run output directory: `{st.session_state.last_run_outdir}`")
    if st.session_state.last_run_outdir and os.path.isdir(st.session_state.last_run_outdir):
        result_files = [f for f in os.listdir(st.session_state.last_run_outdir) if f.endswith((".csv", ".tsv"))]
        if result_files:
            selected_file = st.selectbox("Select a result file to preview", result_files)
            if selected_file:
                df = pd.read_csv(os.path.join(st.session_state.last_run_outdir, selected_file), sep=None, engine='python')
                st.dataframe(df)
        else:
            st.info("No result files (CSV/TSV) found in the output directory.")
    else:
        st.info("Run a job to see results here.")

with tabs[2]:
    if st.session_state.log_text:
        st.code(st.session_state.log_text, language="bash")
    else:
        st.caption("No logs yet. Run a job to see the full log.")

with tabs[3]:
    readme_path = Path(__file__).parent / "README_for_ui.md"
    if readme_path.exists():
        st.markdown(readme_path.read_text(encoding="utf-8"))
    else:
        st.warning("README_for_ui.md not found.")

# --- Footer ---
st.divider()
st.caption(f"AltEx-BE UI | {time.strftime('%Y-%m-%d %H:%M:%S')}")
if st.session_state.last_cmd:
    st.caption("Last command executed:")
    st.code(" ".join(st.session_state.last_cmd), language="bash")
else:
    st.caption("No commands run yet in this session.")