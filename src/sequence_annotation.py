import pandas as pd
import pybedtools

def annotate_sequence_to_bed(bed: pybedtools.BedTool, fasta_path:str,) -> pd.DataFrame:
    """
    Purpose:
        BED形式のデータに指定される遺伝子座位を参照して、塩基配列をFASTAから取得し、bedに塩基配列を追加する
    Parameters:
        bed: pybedtools.Bedtools, BED形式のデータ
        single_exon_df: pd.DataFrame, アノテーションを与えるデータフレーム
    Returns:
        bed_for_df: pd.DataFrame, 配列アノテーションが追加されたデータフレーム(形式はBED)
    """
    bed = bed.copy()
    bed_for_sequence = pybedtools.BedTool(bed)
    # FASTAファイルから配列を取得　s = true　で配列のstrandを考慮し、-の時は相補鎖を出力する
    # もちろん、相補鎖も5'-3'の方向に出力される (今後間違えないように注意する)
    fasta_sequences = bed_for_sequence.sequence(fi=fasta_path, s=True, name=True)
    # 得られた配列ををlistにする(最終的に並べ替えていないBEDに付加するので、keyは必要ない)
    sequences = []
    # ヘッダー行で始まっているときはそれまでに溜めたseqをlistに追加し、seqを初期化
    # それ以外の行は塩基配列としてseqに追加
    with open(fasta_sequences.seqfn) as f:
        seq = ''
        for line in f:
            if line.startswith('>'): #
                if seq: 
                    sequences.append(seq)
                    seq = ''
            else: 
                seq += line.strip() # 改行を除去してseqに追加
        # 以上の処理は次のheaderに到達したときにlistにseqを追加する
        # 最後の配列は次のheaderが存在しないため、for ループを抜けた後にリストに追加する必要がある
        if seq:
            sequences.append(seq)
    # .sequence()メソッドを使用すると、元の構造が保持されないため、bedの別のコピーをdfにする
    bed_for_df = bed.to_dataframe(names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])
    # 取得した配列をデータフレームに追加
    bed_for_df['sequence'] = sequences
    return bed_for_df

def join_sequence_to_single_exon_df(
    single_exon_df: pd.DataFrame,
    bed_for_df: pd.DataFrame,
    acceptor_or_donor: str = 'acceptor' or 'donor', 
    ) -> pd.DataFrame:
    """
    Purpose:
        single_exon_dfにBED形式のデータから、score列取得した塩基配列を追加する
    Parameters:
        single_exon_df: pd.DataFrame, アノテーションを与えるデータフレーム
        bed_for_df: pd.DataFrame, 配列アノテーションが追加されたデータフレーム(形式はBED)
        acceptor_or_donor: str, 'acceptor' or 'donor', どちらの配列を追加するか指定する 
    Returns:
        single_exon_df: pd.DataFrame, 配列アノテーションが追加されたデータフレーム
    """
