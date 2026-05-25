# AltEx-BE: 塩基編集によるエクソンスキッピングを誘導するsgRNA設計ツール

![Build Status](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python Version](https://img.shields.io/badge/python-3.12%2B-blue)
[![PyPI](https://img.shields.io/pypi/v/AltEx-BE)](https://pypi.org/project/AltEx-BE/)

<img src= https://github.com/kinari-labwork/AltEx-BE/raw/refseqid-issue/docs/AltEx-BE_logo.jpg>


- [AltEx-BE: 塩基編集によるエクソンスキッピングを誘導するsgRNA設計ツール](#altex-be-塩基編集によるエクソンスキッピングを誘導するsgrna設計ツール)
  - [重要なお知らせ](#重要なお知らせ)
  - [概要](#概要)
  - [主な特長](#主な特長)
  - [ワークフロー図](#ワークフロー図)
  - [インストール](#インストール)
  - [必要なデータセット](#必要なデータセット)
  - [使い方](#使い方)
    - [グラフィカルユーザーインターフェース (UI)](#グラフィカルユーザーインターフェース-ui)
    - [コマンドラインインターフェース (CLI)](#コマンドラインインターフェース-cli)
      - [1. プリセットエディタを使う場合](#1-プリセットエディタを使う場合)
      - [2. コマンドラインでBase Editorを指定する場合](#2-コマンドラインでbase-editorを指定する場合)
      - [3. Base Editor情報のCSV/TSV/TXTを指定する場合](#3-base-editor情報のcsvtsvtxtを指定する場合)
  - [コマンドラインオプション一覧](#コマンドラインオプション一覧)
  - [AltEx-BEの出力形式](#altex-beの出力形式)
    - [sgRNAの優先順位付け](#sgrnaの優先順位付け)
- [ライセンス](#ライセンス)

## 重要なお知らせ
> [!IMPORTANT]
> **AltEx-BEのWEB版も利用可能です。インストール不要で使えます。**
> 以下のページをご覧ください: https://kinari-labwork.github.io/altex-be-web/

## 概要

**AltEx-BE** は、ベース編集技術を用いて標的エクソンのスキッピングを誘導するためのsgRNA (single guide RNA) を設計するコマンドライン向けバイオインフォマティクスツールです。

代替スプライシングの操作は、がんや神経変性疾患などの理解に不可欠ですが、標的エクソンの同定、ベースエディタに応じたsgRNA設計、オフターゲット評価を手作業で行うのは複雑で時間がかかります。

**AltEx-BE** はこの一連の工程を自動化し、解析の負担を大幅に削減します。refFlat などの注釈から標的候補を抽出し、さまざまなBase Editorに対応したsgRNAを設計し、オフターゲットリスクに基づいて順位付けした候補を出力します。

## 主な特長

- 🧬 **標的エクソンの自動注釈**:
    - refFlatを解析して、Skipped Exon (SE) や Alternative 3'/5' Splice Site (A3SS/A5SS) を自動同定します。

- ⚙️ **多様なBase Editorに対応**:
    - ABE/CBEを問わず、PAM配列と編集ウィンドウを指定するだけで利用できます。
    - 複数のBase Editorを同時に指定して設計可能です。

- 🚀 **エンドツーエンドのワークフロー**:
    - sgRNA設計から出力評価まで一括実行。候補を比較・評価できるレポートを作成します。

## ワークフロー図

**AltEx-BE** の処理フローの概要図です。

<img src = https://github.com/kinari-labwork/AltEx-BE/raw/main/docs/pipeline_explanation.png width="75%">

## インストール

```sh
# 1. bioconda でインストール
conda install -c conda-forge -c bioconda altex-be

# 2. pypi でインストール
pip install AltEx-BE
```

> [!CAUTION]
> AltEx-BEは **Python 3.10~3.12** を必要とします。
> インストールエラーが出る場合は、Pythonのバージョンを確認してください。

## 必要なデータセット
AltEx-BEを実行するには、以下の入力ファイルを準備してください。

- refFlatファイルまたはGTFファイル（対象生物種）
    - refFlatの形式説明: https://genome.bio.fsu.edu/cgi-bin/hgTables?hgsid=235697_cnEhDmy3qVsShD0gwzprkJveBQah&hgta_doSchemaDb=mm39&hgta_doSchemaTable=refFlat
    - refFlatはUCSC Goldenpathから取得できます。例: mm39
      - https://hgdownload.cse.ucsc.edu/goldenpath/mm39/database/
    - GTFも入力として使用可能です（自動でrefFlatへ変換され、出力先に保存されます）

- FASTAファイル（全染色体配列）
    - UCSC Goldenpathから取得可能です。
    - `.fa`に全染色体が含まれていない場合、処理が失敗します。

- (任意) 遺伝子シンボルまたはRefSeq IDの一覧ファイル（CSV/TXT/TSV）
    - 複数遺伝子に対して一括設計したい場合に `--gene-file` で指定できます。
    - 1列のみ（ヘッダー不要）で、遺伝子シンボルまたはRefSeq IDを記載してください。

> [!NOTE]
> **遺伝子シンボルとRefSeq IDの入力について**
> - 遺伝子シンボル（例: MYGENE）を指定した場合、その遺伝子の既知転写産物全体を解析し、代替スプライシングイベントを探索します。
> - RefSeq ID（例: NM_0012345）を指定した場合、自動的に対応遺伝子を特定し、同様に包括的な解析を行います。

## 使い方

AltEx-BEはGUIまたはCLIから利用できます。

### グラフィカルユーザーインターフェース (UI)

StreamlitベースのWeb UIを起動して、コマンドラインを使わずに操作できます。

**起動方法:**
```sh
altex-be --ui
```

**UIでできること:**
- 📥 入力ファイルの選択（FASTA、注釈）
- ⚙️ オプション設定
- ▶️ パイプライン実行
- 📝 ログのリアルタイム表示
- 📊 出力ファイルのプレビュー

### コマンドラインインターフェース (CLI)

AltEx-BEは `altex-be` コマンドで実行します。例:

#### 1. プリセットエディタを使う場合

> [!NOTE]
> **プリセットBase Editor:**
>
> | base_editor_name   | pam_sequence | editing_window_start | editing_window_end | base_editor_type |
> |:------------------|:-------------|:---------------------|:-------------------|:-----------------|
> | target_aid_ngg     | NGG          | 17                   | 19                 | cbe              |
> | be4max_ngg         | NGG          | 12                   | 17                 | cbe              |
> | abe8e_ngg          | NGG          | 12                   | 17                 | abe              |
> | target_aid_ng      | NG           | 17                   | 19                 | cbe              |
> | be4max_ng          | NG           | 12                   | 17                 | cbe              |
> | abe8e_ng           | NG           | 12                   | 17                 | abe              |

```sh
altex-be \
    --refflat-path /path/to/your/refFlat.txt \
    --fasta-path /path/to/your/genome.fa \
    --output-dir /path/to/output_directory \
    --gene-symbols MYGENE \
    --assembly-name hg38
```

#### 2. コマンドラインでBase Editorを指定する場合

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
> `--be-start` と `--be-end` は編集ウィンドウの範囲です。
> PAM直後の塩基を1番としてカウントします。

#### 3. Base Editor情報のCSV/TSV/TXTを指定する場合

> [!CAUTION]
> 入力ファイルは以下の列を含めてください: `base_editor_name`, `pam_sequence`, `editing_window_start`, `editing_window_end`, `base_editor_type`.

```sh
altex-be \
    --refflat-path /path/to/your/refFlat.txt \
    --fasta-path /path/to/your/genome.fa \
    --output-dir /path/to/output_directory \
    --gene-symbols MYGENE \
    --assembly-name hg38 \
    --be-files /path/to/your/base_editor_info.csv
```

## コマンドラインオプション一覧

| Short Option | Long Option | Argument | Explanation |
| :--- | :--- | :--- | :--- |
| -h | --help | | ヘルプを表示して終了します。 |
| -v | --version | | AltEx-BEのバージョンを表示します。 |
| | --ui | | Streamlit UIを起動します。 |
| -r | --refflat-path | FILE | (必須: -r または -g) refFlatファイルのパス |
| -g | --gtf-path | FILE | (必須: -r または -g) GTFファイルのパス |
| -f | --fasta-path | FILE | (必須) FASTAファイルのパス |
| -o | --output-dir | DIR | (必須) 出力ディレクトリ |
| | --gene-symbols| SYMBOL [SYMBOL ...] | 対象遺伝子シンボルのリスト |
| | --refseq-ids | ID [ID ...] | 対象RefSeq IDのリスト |
| | --gene-file | FILE | 遺伝子シンボル/RefSeq ID一覧ファイル |
| -a | --assembly-name| ASSEMBLY | (必須) ゲノムアセンブリ名（例: hg38, mm39） |
| -n | --be-name | NAME | Base Editor名 |
| -p | --be-pam | SEQUENCE | PAM配列 |
| -s | --be-start | INTEGER | 編集ウィンドウ開始位置 (1-indexed) |
| -e | --be-end | INTEGER | 編集ウィンドウ終了位置 (1-indexed) |
| -t | --be-type | TYPE | Base Editorタイプ（ABE/CBE） |
| | --be-files | FILE | Base Editor情報ファイル（CSV/TXT） |

## AltEx-BEの出力形式

`altex-be` は `--output-dir` で指定したディレクトリに2つの出力ファイルを生成します。

- sgRNAサマリーテーブル (.csv)
    - 設計されたsgRNAの情報と順位が含まれます。
<img src=https://github.com/kinari-labwork/AltEx-BE/raw/main/docs/output_csv_example.png width = "100%">

**各列の意味:**

|column name|meaning|remark|
|:-----------|:-------------------|--|
|geneName|対象遺伝子シンボル| |
|chrom|対象遺伝子の染色体| |
|strand|対象遺伝子のストランド| |
|exonstart, exonend, exonlength|標的エクソンの情報| |
|coding|タンパク質コードか非コードか| |
|frame|エクソン長のmod3|0 = in-frame / 1,2 = out-frame |
|exon_position|標的エクソンの位置|"first" / "internal" / "last"|
|uuid|sgRNAの一意ID|実行ごとに変わる|
|exon_intron_boundary+-25bp_sequence|SA/SD周辺配列| |
|sgrna_sequence|sgRNA配列|TはUに置換しない|
|sgrna_target_pos_in_seq|sgRNA中の編集塩基位置|相対位置|
|sgrna_overlap_between_cds_and_editing_window|編集ウィンドウとCDSの重複塩基数| |
|sgrna_unintended_edited_base_count|CDS内の編集候補塩基数| |
|sgrna_start/end_in_genome|sgRNAのゲノム座標| |
|site type|標的スプライシング部位|acceptor / donor|
|sgrna_strand|sgRNAのストランド| |
|base_editor_name/pam_sequence/window_start or end / base editor type|Base Editor情報| |
|crispr_direct_url|CRISPR directへのリンク| |
|pam+20bp exact match|PAM+20bpの完全一致数| |
|pam+12bp exact match|PAM+12bpの完全一致数| |
|sgrna_priority|sgRNA順位|オフターゲットとGCで評価|

### sgRNAの優先順位付け

AltEx-BEは、各スプライス部位に対してオフターゲット特異性が高いsgRNAを優先します。まず、ゲノム中のPAM+20bpの完全一致数が少ないものを優先し、同点の場合はGC含量（40-60%が目安）で順位付けします。それでも同点の場合は、PAM+12bp一致数やCDS内編集候補塩基数をタイブレークに使用します。

- UCSCカスタムトラック用BED (.bed)
    - UCSCカスタムトラックとして利用できます: https://genome.ucsc.edu/cgi-bin/hgCustom
<img src = https://github.com/kinari-labwork/AltEx-BE/raw/main/docs/examle_of_custom_track.png width = "75%">
    - 赤/青のボックスはsgRNA配列です（赤=ABE, 青=CBE）。
    - BEDのscore列は20bp+PAMのオフターゲット数を示します。
    - 使用時はUCSC側で正しいアセンブリ名を選択してください。

# ライセンス
- 詳細は [LICENSE.md](LICENSE.md) を参照してください。
