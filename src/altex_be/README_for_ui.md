# 🧬 AltEx-BE (UI)

AltEx-BE UI is a **local graphical interface** for the AltEx-BE command-line tool.
It helps you configure inputs and options, run the analysis, and preview results
without writing long command-line arguments.

---

## What this UI does

- 📥 Select input files (FASTA, transcript annotation)
- ⚙️ Configure AltEx-BE command-line options
- ▶ Run AltEx-BE locally
- 📝 Monitor execution logs in real time
- 📊 Preview output tables (CSV / TSV)
- 📁 Browse generated output files

> This UI is intended for **local use** after installing AltEx-BE via `conda` or `pip`.

---

## Typical workflow

1. Select **FASTA** and **annotation file** in the sidebar  
2. Choose an **output directory**
3. Configure options (mode, threads, base editors, etc.)
4. Click **Run AltEx-BE**
5. Monitor progress in the log panel
6. Preview results in the **Results** tab

---

## Notes

- This UI is a wrapper around the `altex-be` CLI  
- All computations are performed locally
- For detailed algorithm descriptions and full CLI usage, see the full README:

👉 **Full documentation:**  
https://github.com/kinari-labwork/AltEx-BE

---

## Citation / License

AltEx-BE is distributed under the MIT License.  
Please cite the associated publication when using this tool in academic work.
