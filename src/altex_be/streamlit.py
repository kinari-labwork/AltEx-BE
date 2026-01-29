import streamlit as st
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import os

# streamlit 1.53.1

apptitle = "Altex BE"
st.set_page_config(page_title=apptitle, layout="wide")

df = pd.DataFrame({
    'test1' : [1, 2, 3, 4],
    'test2' : [10, 20, 30, 40]
})
st.write(df)

def select_folder():
   root = tk.Tk()
   root.withdraw()
   folder_path = filedialog.askdirectory(
    master=root)
   root.destroy()
   return folder_path

def select_file():
   root = tk.Tk()
   root.withdraw()
   file_path = filedialog.askopenfilename(
    master=root,
    filetypes=[("Fasta files", "*.fa"), ("All files", "*.*")])
   root.destroy()
   return file_path

def is_writable_dir(path):
    return os.path.isdir(path) and os.access(path, os.W_OK)


if st.button('Select Output Directory'):
    folder = select_folder()
    if folder:
        st.write(f'Selected Output Directory: {folder}')
    else:
        st.write('No Output Directory selected.')

if st.button('Select Input Fasta File'):
    file = select_file()
    if file:
        st.write(f'Selected Input Fasta File: {file}')
    else:
        st.write('No Input Fasta File selected.')