import pandas as pd
import pybedtools



def annotate_sequence_to_bed(bed: pd.DataFrame, fasta_path: str) -> pd.DataFrame:
    """
    Purpose:
        BED形式のデータに指定される遺伝子座位を参照して、塩基配列をFASTAから取得し、bedに塩基配列を追加する
    Parameters:
        bed: pd.DataFrame, BED形式のデータ(内部でpybedtoolsを使用するため、pybedtools.BedToolに変換可能な形式である必要がある,
            例: ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']の列を持つ)
        fasta_path: str, FASTAファイルのパス
    Returns:
        bed_for_df: pd.DataFrame, 配列アノテーションが追加されたデータフレーム(形式はBED)
    """
    bed_for_df = bed
    bed_for_sequence = pybedtools.BedTool.from_dataframe(bed)
    # FASTAファイルから配列を取得　s = true　で配列のstrandを考慮し、-の時は相補鎖を出力する
    # もちろん、相補鎖も5'-3'の方向に出力される (今後間違えないように注意する)
    fasta_sequences = bed_for_sequence.sequence(fi=fasta_path, s=True, name=True)
    # 得られた配列ををlistにする(最終的に並べ替えていないBEDに付加するので、keyは必要ない)
    sequences = []
    # ヘッダー行で始まっているときはそれまでに溜めたseqをlistに追加し、seqを初期化
    # それ以外の行は塩基配列としてseqに追加
    with open(fasta_sequences.seqfn) as f:
        seq = []
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append("".join(seq))
                    seq = []
            else:
                seq.append(line.strip())  # 改行を除去してseqに追加
        # 以上の処理は次のheaderに到達したときにlistにseqを追加する
        # 最後の配列は次のheaderが存在しないため、for ループを抜けた後にリストに追加する必要がある
        if seq:
            sequences.append("".join(seq))
    # .sequence()メソッドを使用すると、元の構造が保持されないため、bedの別のコピーをdfにする
    bed_for_df.columns = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    # 取得した配列をデータフレームに追加
    bed_for_df["sequence"] = sequences
    return bed_for_df


def join_sequence_to_single_exon_df(
    single_exon_df: pd.DataFrame,
    acceptor_bed_with_sequences: pd.DataFrame,
    donor_bed_with_sequences: pd.DataFrame,
) -> pd.DataFrame:
    """
    Purpose:
        single_exon_dfにbed_for_dfで取得した塩基配列をleft joinして、acceptor_sequenceまたはdonor_sequence列を追加する
    Parameters:
        single_exon_df: pd.DataFrame, アノテーションを与えるデータフレーム
        acceptor_bed_with_sequences: pd.DataFrame, acceptorのBED形式のデータフレームに配列アノテーションが追加されたもの
        donor_bed_with_sequences: pd.DataFrame, donorのBED形式のデータフレーム
    Returns:
        single_exon_df: pd.DataFrame, 配列アノテーションが追加されたデータフレーム
    """

    acceptor_bed_with_sequences = acceptor_bed_with_sequences.rename(
        columns = {"name":"uuid","chromStart": "chromStart_acceptor", "chromEnd": "chromEnd_acceptor"})
    
    donor_bed_with_sequences = donor_bed_with_sequences.rename(
        columns = {"name":"uuid","chromStart": "chromStart_donor", "chromEnd": "chromEnd_donor"})

    for label, bed in [
        ("acceptor", acceptor_bed_with_sequences),
        ("donor", donor_bed_with_sequences),
    ]:
        single_exon_df = single_exon_df.merge(
            bed[["uuid", f"chromStart_{label}", f"chromEnd_{label}", "sequence"]], on="uuid", how="left"
        ).rename(columns={"sequence": f"{label}_exon_intron_boundary_±25bp_sequence"})
    return single_exon_df

def annotate_sequence_to_splice_sites(
    single_exon_df: pd.DataFrame,
    splice_acceptor_single_exon_df: pd.DataFrame,
    splice_donor_single_exon_df: pd.DataFrame,
    fasta_path: str
) -> pd.DataFrame:
    """
    このモジュールの操作をまとめて実行するためのラッパー関数
    """
    acceptor_bed_with_sequences = annotate_sequence_to_bed(
        splice_acceptor_single_exon_df, fasta_path
    )
    donor_bed_with_sequences = annotate_sequence_to_bed(
        splice_donor_single_exon_df, fasta_path
    )
    single_exon_df = join_sequence_to_single_exon_df(single_exon_df, acceptor_bed_with_sequences, donor_bed_with_sequences)
    return single_exon_df
