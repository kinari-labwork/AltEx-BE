import pandas as pd
from altex_aid.sgrna_designer import (
    design_sgrna_for_target_exon_df,
    organize_target_exon_df_with_grna_sequence,
    convert_sgrna_start_end_position_to_position_in_chromosome,
)

target_exons_df_with_acceptor_and_donor_sequence = pd.read_pickle("data/target_exons_with_acceptor_and_donor_sequence.pkl")


# とりあえず target-AIDの条件で実行してみる
target_exons_df_with_grna_info = design_sgrna_for_target_exon_df(
    target_exon_df=target_exons_df_with_acceptor_and_donor_sequence,
    pam_sequence="NGG",
    editing_window_start_in_grna=17,
    editing_window_end_in_grna=19
)
target_exons_df_with_grna_info = organize_target_exon_df_with_grna_sequence(target_exons_df_with_grna_info)

target_exons_df_with_grna_info = convert_sgrna_start_end_position_to_position_in_chromosome(target_exons_df_with_grna_info)

target_exons_df_with_grna_info.to_pickle("data/target_exons_with_grna_sequence.pkl")