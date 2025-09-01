# AltexBE: Alternate Exon Skipping by Base Editing

![Build Status](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python Version](https://img.shields.io/badge/python-3.12%2B-blue)

- [AltexBE: Alternate Exon Skipping by Base Editing](#altexbe-alternate-exon-skipping-by-base-editing)
  - [Overview](#overview)
  - [Key Features](#key-features)
  - [Workflow Diagram](#workflow-diagram)
  - [Installation](#installation)
  - [Required dataset](#required-dataset)
  - [Format of Altex-BE output](#format-of-altex-be-output)
  - [Usage](#usage)
    - [List of command line options](#list-of-command-line-options)
      - [1. Input Base Editor Information in the Command Line:](#1-input-base-editor-information-in-the-command-line)
      - [2. Input a CSV/TSV/TXT File Containing Information about Your Base Editors:](#2-input-a-csvtsvtxt-file-containing-information-about-your-base-editors)
      - [3. Using a Preset Editor:](#3-using-a-preset-editor)
  - [License](#license)

## Overview

**AltexBE** is a command-line bioinformatics tool that designs sgRNAs (single guide RNAs) to induce targeted exon skipping using Base Editing technology.

Manipulating alternative splicing is key to understanding diseases like cancer and neurodegenerative disorders, but designing the right tools for the job is a major bottleneck. The manual process of identifying targetable exons, designing sgRNAs for specific base editors, and assessing off-target risks is complex, tedious, and slows down critical research.

**AltEx-BE** is a powerful command-line tool built to automate this entire workflow. It intelligently parses transcript data to find the best exon targets, designs candidates for a multitude of base editors, and evaluates their off-target risk to provide a ranked list of high-confidence sgRNAs.

By transforming a complex, multi-step design process into a single command, AltEx-BE bridges the gap between your scientific question and a successful wet lab experiment, significantly accelerating research into splicing-related diseases and therapies.


## Key Features

- 🧬 **Automated Target Exon Annotation**: 
    - Automatically parses transcript structures from refFlat files to identify and classify potential targets for exon skipping. This includes Skipped Exons (SE) and exons with Alternative 3'/5' Splice Sites (A3SS/A5SS), eliminating the need for tedious manual searches.

- ⚙️ **Universal Base Editor Compatibility**: 
    - Supports virtually any ABE or CBE. You can use built-in presets or define any custom editor by specifying its PAM sequence and editing window, allowing immediate use of the latest editors from new publications. 
    - AltEx-BE can design sgRNAs for multiple Base-Editors in one run

- 🚀 **Streamlined End-to-End Workflow**:
    - Seamlessly moves from data input to candidate selection. The design command generates sgRNAs, while the visualize command creates comprehensive reports to help you evaluate and rank the best candidates for your experiment.

## Workflow Diagram

Here is a simplified diagram illustrating the workflow of `AltexBE`:

![pipeline_explanation](docs/pipeline_explanation.png)

## Installation

To get started with AltexBE, clone the repository and install the required dependencies.

```sh
# 1. Clone the repository
git clone https://github.com/kinari-labwork/AltEx-BE
cd AltEx-BE
pip install -e .
```
## Required dataset
To use AltexBE, you should prepare 2 input files in your computer
- refFlat file of your interest species   
    - refflat file contains Refseq infomations: explanation of refFlat format is [here](https://genome.bio.fsu.edu/cgi-bin/hgTables?hgsid=235697_cnEhDmy3qVsShD0gwzprkJveBQah&hgta_doSchemaDb=mm39&hgta_doSchemaTable=refFlat)   
    - you can download refflat files from  UCSC goldenpath: refflat files of mm39 is [here](https://hgdownload.cse.ucsc.edu/goldenpath/mm39/database/)
- Fasta files contain all chromosome sequence of your interest species
    - you can download Fasta file also from UCSC goldenpath
    - please comfirm your .fa files contain all of chromosome. if not, AltEx-BE process will fail

## Format of Altex-BE output
`altex-be` makes 2 output files in `Path/To/YourOutput/` directory which you specified in `--output-dir` command
- Summary sgRNA table (.csv)
    - this table contain imformation of sgRNAs designed by AltEx-BE
- BED file for UCSC custom track (.bed)
    - this bed file can use as a UCSC custom tracks, you can input that bed file into [this webpage](https://genome.ucsc.edu/cgi-bin/hgCustom)
    - when you assign bed file, you should choose correct assembly name in above website
    - score columns in bed file means offtarget count of 20bp+PAM

## Usage

AltexBE is operated via the `altex-be` command.

### List of command line options

| Short Option | Long Option | Argument | Explanation |
| :--- | :--- | :--- | :--- |
| -h | --help | | Show the help message and exit. |
| -v | --version | | Show the version of Altex BE. |
| -r | --refflat-path | FILE | (Required) Path to the refFlat file. |
| -f | --fasta-path | FILE | (Required) Path to the FASTA file. |
| -o | --output-dir | DIR | (Required) Directory for the output files. |
| | --gene-symbols| SYMBOL [SYMBOL ...] | A space-separated list of gene symbols of interest. |
| | --refseq-ids | ID [ID ...] | A space-separated list of RefSeq IDs of interest. |
| -a | --assembly-name| ASSEMBLY | (Required) The name of the genome assembly to use (e.g., hg38, mm39). |
| -n | --be-name | NAME | The name of the base editor to use. |
| -p | --be-pam | SEQUENCE | The PAM sequence for the base editor. |
| -s | --be-start | INTEGER | The start of the editing window for the base editor (1-indexed from the base next to the PAM). |
| -e | --be-end | INTEGER | The end of the editing window for the base editor (1-indexed from the base next to the PAM). |
| -t | --be-type | TYPE | The type of base editor (ABE or CBE). |
| | --be-preset | PRESET | Use a preset base editor (target-AID, BE4max, or ABE8e). |
| | --be-files | FILE | Path to a CSV or TXT file containing information about one or more base editors. |

#### 1. Input Base Editor Information in the Command Line:

```sh
altex-be \
    --refflat-path /path/to/your/refFlat.txt \
    --fasta-path /path/to/your/genome.fa \
    --output-dir /path/to/output_directory \
    --gene-symbols MYGENE \
    --assembly-name hg38 \
    --be-name target-aid \
    --be-type cbe \
    --be-pam NGG \
    --be-start 17 \
    --be-end 19
```

> [!CAUTION]
> `--be-start` and `--be-end` specify the editing window of your base editor. The location of the editing window is counted from the base next to the PAM (1-indexed).

#### 2. Input a CSV/TSV/TXT File Containing Information about Your Base Editors:

You can provide a file containing the information for one or more base editors. This is useful when you want to design sgRNAs for multiple editors at once.

The input file should have the following columns: `base_editor_name`, `pam_sequence`, `editing_window_start`, `editing_window_end`, `base_editor_type`.

```sh
altex-be \
    --refflat-path /path/to/your/refFlat.txt \
    --fasta-path /path/to/your/genome.fa \
    --output-dir /path/to/output_directory \
    --gene-symbols MYGENE \
    --assembly-name hg38 \
    --be-files /path/to/your/base_editor_info.csv
```

#### 3. Using a Preset Editor:

You can use a pre-configured base editor with the `--be-preset` flag.

> [!NOTE]
> **Preset Base Editors:**
>
> | base_editor_name | pam_sequence | editing_window_start | editing_window_end | base_editor_type |
> |:-----------------|:-------------|:---------------------|:-------------------|:-----------------|
> | target-AID       | NGG          | 17                   | 19                 | cbe              |
> | BE4max           | NGG          | 12                   | 17                 | cbe              |
> | ABE8e            | NGG          | 12                   | 17                 | abe              |

```sh
altex-be \
    --refflat-path /path/to/your/refFlat.txt \
    --fasta-path /path/to/your/genome.fa \
    --output-dir /path/to/output_directory \
    --gene-symbols MYGENE \
    --assembly-name hg38 \
    --be-preset ABE8e
```

## License

This project is distributed under the MIT License. See the `LICENSE` file for more information.
