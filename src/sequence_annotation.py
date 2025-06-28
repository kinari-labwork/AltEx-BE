import pandas as pd
from typing import Literal
import pybedtools


# 今はまだFasta pathがハードコーディングになっているので、再利用性を考えて後でargpersで指定できるようにする
def annotate_sequence_to_bed(bed: pd.DataFrame, fasta_path:str) -> pd.DataFrame:
    """
    Purpose:
        BED形式のデータに指定される遺伝子座位を参照して、塩基配列をFASTAから取得し、bedに塩基配列を追加する
    Parameters:
        bed: pd.DataFrame, BED形式のデータ(内部でpybetdtoolsを使用するため、pybedtools.BedToolに変換可能な形式である必要がある,
            例: ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']の列を持つ)
        fasta_path: str, FASTAファイルのパス
    Returns:
        bed_for_df: pd.DataFrame, 配列アノテーションが追加されたデータフレーム(形式はBED)
    """
    bed = bed.copy()
    bed_for_df = bed.copy()
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
            if line.startswith('>'):
                if seq: 
                    sequences.append(''.join(seq))
                    seq = []
            else: 
                seq.append(line.strip())# 改行を除去してseqに追加
        # 以上の処理は次のheaderに到達したときにlistにseqを追加する
        # 最後の配列は次のheaderが存在しないため、for ループを抜けた後にリストに追加する必要がある
        if seq:
            sequences.append(''.join(seq))
    # .sequence()メソッドを使用すると、元の構造が保持されないため、bedの別のコピーをdfにする
    bed_for_df.columns = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
    # 取得した配列をデータフレームに追加
    bed_for_df['sequence'] = sequences
    return bed_for_df

def join_sequence_to_single_exon_df(
    single_exon_df: pd.DataFrame,
    acceptor_or_donor_bed_with_sequences: pd.DataFrame,
    acceptor_or_donor: Literal['acceptor', 'donor']
    ) -> pd.DataFrame:
    """
    Purpose:
        single_exon_dfにbed_for_dfで取得した塩基配列をleft joinして、acceptor_sequenceまたはdonor_sequence列を追加する
    Parameters:
        single_exon_df: pd.DataFrame, アノテーションを与えるデータフレーム
        bed_for_df: pd.DataFrame, 配列アノテーションが追加されたデータフレーム(形式はBED)
        acceptor_or_donor: str, 'acceptor' or 'donor', どちらの配列を追加するか指定する
    Returns:
        single_exon_df: pd.DataFrame, 配列アノテーションが追加されたデータフレーム
    """
    # acceptor_or_donorに応じて、single_exon_dfに加える列名を決定
    if acceptor_or_donor == 'acceptor':
        column_name = 'acceptor_sequence'
    elif acceptor_or_donor == 'donor':
        column_name = 'donor_sequence'
    else:
        raise ValueError("acceptor_or_donor must be 'acceptor' or 'donor'")
    # single_exon_dfのscore列をキーとして、bed_for_dfのsequence列を[column_name]列に追加
    single_exon_df = single_exon_df.merge(
        acceptor_or_donor_bed_with_sequences[['score','sequence']],
        on = 'score',
        how ='left',
        ).rename(columns={'sequence': column_name})
    return single_exon_df

