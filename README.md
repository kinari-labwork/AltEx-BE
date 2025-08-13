# AltexBE: Alternate Exon Skipping by Base Editing

![Build Status](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python Version](https://img.shields.io/badge/python-3.9%2B-blue)

## Overview

**AltexBE** is a command-line bioinformatics tool that designs sgRNAs (single guide RNAs) to induce targeted exon skipping using Base Editing technology.

Manipulating alternative splicing is key to understanding diseases like cancer and neurodegenerative disorders, but designing the right tools for the job is a major bottleneck. The manual process of identifying targetable exons, designing sgRNAs for specific base editors, and assessing off-target risks is complex, tedious, and slows down critical research.

**AltEx-BE** is a powerful command-line tool built to automate this entire workflow. It intelligently parses transcript data to find the best exon targets, designs candidates for a multitude of base editors, and evaluates their off-target risk to provide a ranked list of high-confidence sgRNAs.

By transforming a complex, multi-step design process into a single command, AltEx-BE bridges the gap between your scientific question and a successful wet lab experiment, significantly accelerating research into splicing-related diseases and therapies.


## Key Features

- 🧬 Automated Target Exon Annotation: 
    - Automatically parses transcript structures from refFlat files to identify and classify potential targets for exon skipping. This includes Skipped Exons (SE) and exons with Alternative 3'/5' Splice Sites (A3SS/A5SS), eliminating the need for tedious manual searches.

- ⚙️ Universal Base Editor Compatibility: 
    - Supports virtually any ABE or CBE. You can use built-in presets or define any custom editor by specifying its PAM sequence and editing window, allowing immediate use of the latest editors from new publications. 
    - AltEx-BE can design sgRNAs for multiple Base-Editors in one run

- 🚀 Streamlined End-to-End Workflow:
    - Seamlessly moves from data input to candidate selection. The design command generates sgRNAs, while the visualize command creates comprehensive reports to help you evaluate and rank the best candidates for your experiment.

## Workflow Diagram

Here is a simplified diagram illustrating the workflow of `AltexBE`:

```
[論文のfig1に相当する部分を入れる予定です]
(e.g., refFlat + FASTA -> AltexBE `design` -> sgRNA candidates (pickle) -> AltexBE `visualize` -> )
```

## Installation

To get started with AltexBE, clone the repository and install the required dependencies.

```sh
# 1. Clone the repository
git clone https://github.com/kinari-labwork/AltEx-BE
pip install -e .

pip install altex-be
```

## Usage

AltexBE is operated via two main subcommands: `design` and `visualize`.

---

## 1. `design` Command

This command designs sgRNAs based on your input files and saves the results as a csv file.  
There are 3 ways to input your interest base editing tools

---

### **Input base editor infomation in command line:**

Here is a basic command to design sgRNAs for the gene *MYGENE*.

```sh
altex-be design \
    --refflat-path /path/to/your/refFlat.txt \
    --fasta-path /path/to/your/genome.fa \
    --output-dir /path/to/output_directory \
    --gene-symbols MYGENE \
    --base-editor-name target-aid \
    --base-editor-type cbe \
    --base-editor-window-start 17 \
    --base-editor-window-end 19 \
```

> [!CAUTION]
> `--base-editor-window-start/end` means editing window of your base editing tool.   
> **Location of editing window Start and END is counted from the base next to PAM (1-index)**

---

### **Input csv, txt, tsv containing infomation of your base-editors**

> [!IMPORTANT]
> AltEx-BE can design sgRNA for multiple Base-Editors
> If you want to do that, you should use following code
> Input `base_editor_info.csv` should follow below format
> |base_editor_name|pam_seqence|editing_window_start|editing_window_end|base_editor_type|
> |---|---|---|---|---|
> |your BE name|Any Sequence|1<int<20|1<int<20|cbe or abe|

```sh
altex-be design \
    --refflat-path /path/to/your/refFlat.txt \
    --fasta-path /path/to/your/genome.fa \
    --output-dir /path/to/output_directory \
    --gene-symbols MYGENE \
    --base-editor-file /path/to/your/base_editor_info.csv
```

---

### **Using a Preset Editor:**

You can specify a pre-configured base editor using the `--editor-preset` flag.
> [!NOTE]
> Preset Base Editors:
> |base_editor_name|pam_seqence|editing_window_start|editing_window_end|base_editor_type|
> |---|---|---|---|---|
> |target-AID|CBE|NGG|17|19|
> |BE4max|CBE|NGG|12|17|
> |ABE8e|ABE|NGG|12|17|

```sh
altex-be design \
    --refflat-path /path/to/your/refFlat.txt \
    --fasta-path /path/to/your/genome.fa \
    --output-dir /path/to/output_directory \
    --gene-symbols MYGENE \
    --editor-preset ABE8e
```

---

### 2. `visualize` Command 仮で置いてます

This command takes the pickle file generated by `design` and creates a visual report.

```sh
python altexbe/main.py visualize \
    --input-pickle /path/to/output_directory/results.pickle \
    --output-report /path/to/output_directory/report
```

## Notes & Warnings

> [!NOTE]
> **refFlat File Format**
> The input `refFlat.txt` file must adhere to the standard format. You can obtain this file for your genome of interest from the UCSC Genome Browser. 
> **reference genome Fasta File Format**
> The imformation of reference genome `your_interest_spiecies.fa` is also required. You can also obtain genome Fasta File from UCSC genome Browser.

> [!WARNING]
> The input `your_interest_spiecies.fa` should contain imformation of all chromosome. You can obtain combined Fasta from UCSC genome browser. 


## License

This project is distributed under the MIT License. See the `LICENSE` file for more information.