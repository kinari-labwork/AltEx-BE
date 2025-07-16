import pandas as pd
from altex_aid.sgrna_designer import (
    SgrnaInfo,
    get_reversed_complement,
    reverse_pam_sequence,
    design_sgrna,
    design_sgrna_for_target_exon_df,
    extract_sgrna_features,
    organize_target_exon_df_with_grna_sequence,
    modify_sgrna_start_end_position_to_position_in_chromosome
)   

pd.set_option('display.max_columns', None)  # 全てのカラムを表示するための設定

def test_get_reversed_complement():
    input_sequence = "ATGCATGC"
    expected_output = "GCAUGCAU" # RNAにした逆相補鎖
    output_sequence = get_reversed_complement(input_sequence)
    assert output_sequence == expected_output

def test_reverse_pam_sequence():
    input_sequence = "NGG"
    expected_output = "[Cc][Cc][ATGCatgc]" 
    output_sequence = reverse_pam_sequence(input_sequence)
    assert output_sequence == expected_output

def test_design_sgrna():
    editing_sequence = "NNNCCCCCNNNNNNNNNNNNNNNAGGGNNNNNNNNNNNNNNNNNNNNNNN"
    # これは+ strandのエキソンで、SAの配列。
    pam_sequence = "NGG"
    editing_window_start_in_grna = 17
    editing_window_end_in_grna = 19
    target_g_pos_in_sequence = 24 # acceptorなら24番目のG, donorなら25番目のGが編集ターゲット
    cds_boundary = 25
    site_type = "acceptor"

    expected_output = [
        SgrnaInfo(
            target_sequence="CCNNNNNNNNNNNNNNNAGG",
            actual_sequence="CCUNNNNNNNNNNNNNNNGG",
            start_in_sequence = 6,
            end_in_sequence=26,
            target_pos_in_sgrna=19,
            overlap_between_cds_and_editing_window=0,
            possible_unintended_edited_base_count= 0
        ),
        SgrnaInfo(
            target_sequence="CNNNNNNNNNNNNNNNAGGG",
            actual_sequence="CCCUNNNNNNNNNNNNNNNG",
            start_in_sequence = 7,
            end_in_sequence=27,
            target_pos_in_sgrna=18,
            overlap_between_cds_and_editing_window=1,
            possible_unintended_edited_base_count= 1
        ),
        SgrnaInfo(
            target_sequence="NNNNNNNNNNNNNNNAGGGN",
            actual_sequence="NCCCUNNNNNNNNNNNNNNN",
            start_in_sequence = 8,
            end_in_sequence=28,
            target_pos_in_sgrna=17,
            overlap_between_cds_and_editing_window=2,
            possible_unintended_edited_base_count= 2
        ),
]
    output_data = design_sgrna(
        editing_sequence=editing_sequence,
        pam_sequence=pam_sequence,
        editing_window_start_in_grna=editing_window_start_in_grna,
        editing_window_end_in_grna=editing_window_end_in_grna,
        target_g_pos_in_sequence=target_g_pos_in_sequence,
        cds_boundary=cds_boundary,
        site_type=site_type
    )
    print(output_data)
    assert output_data == expected_output

def test_design_sgrna_for_target_exon_df():
    target_exon_df = pd.DataFrame({
        "chrom" : ["chr1", "chr2"],
        "exonStarts": [100, 300],
        "exonEnds": [200, 400],
        "name":["UUID1", "UUID2"],
        "score":[0, 0],
        "strand":["+","-"],
        "exontype":["skipped", "a5ss-long"],
        "exon_position": ["internal", "internal"],
        "chromStart_acceptor":[75, 275],
        "chromEnd_acceptor":[125, 325],
        "acceptor_sequence":["NNNNCCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNN","NNNNCCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNN"],
        "chromStart_donor":[175, 375],
        "chromEnd_donor":[225, 425],
        "donor_sequence":["NNNNNCCCNNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNNNNNNNNNNNN","NNNNNCCCNNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNNNNNNNNNNNN"]
    })
    expected_output = pd.DataFrame({
        "chrom": ["chr1", "chr2"],
        "exonStarts": [100, 300],
        "exonEnds": [200, 400],
        "name": ["UUID1", "UUID2"],
        "score": [0, 0],
        "strand": ["+", "-"],
        "exontype": ["skipped", "a5ss-long"],
        "exon_position": ["internal", "internal"],
        "chromStart_acceptor": [75, 275],
        "chromEnd_acceptor": [125, 325],
        "acceptor_sequence":["NNNNCCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNN","NNNNCCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNN"],
        "chromStart_donor": [175, 375],
        "chromEnd_donor": [225, 425],
        "donor_sequence":["NNNNNCCCNNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNNNNNNNNNNNN","NNNNNCCCNNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNNNNNNNNNNNN"],
        "grna_acceptor": [
            [
            SgrnaInfo(
                target_sequence="NNNNNNNNNNNNNNNNAGNN",
                actual_sequence="NNCUNNNNNNNNNNNNNNNN",
                start_in_sequence=7,
                end_in_sequence=27,
                target_pos_in_sgrna=18,
                overlap_between_cds_and_editing_window=1,
                possible_unintended_edited_base_count=0
            )
            ],
            []
        ],
        "grna_donor": [
            [
            SgrnaInfo(
                target_sequence="NNNNNNNNNNNNNNNNNGTN",
                actual_sequence="NACNNNNNNNNNNNNNNNNN",
                start_in_sequence=8,
                end_in_sequence=28,
                target_pos_in_sgrna=18,
                overlap_between_cds_and_editing_window=1,
                possible_unintended_edited_base_count=0
            )
            ],
            [
            SgrnaInfo(
                target_sequence="NNNNNNNNNNNNNNNNNGTN",
                actual_sequence="NACNNNNNNNNNNNNNNNNN",
                start_in_sequence=8,
                end_in_sequence=28,
                target_pos_in_sgrna=18,
                overlap_between_cds_and_editing_window=1,
                possible_unintended_edited_base_count=0
            )
            ]
        ]
        })
    output_data = design_sgrna_for_target_exon_df(
        target_exon_df =target_exon_df,
        pam_sequence="NGG",
        editing_window_start_in_grna=17,
        editing_window_end_in_grna=19,
        )
    print(output_data["grna_acceptor"].tolist())
    print(output_data["grna_donor"].tolist())

    # pd.testing.assert_frame_equalの問題で、SgrnaInfoオブジェクトのようなカスタムクラスを直接比較するとエラーが出るため、リストに変換して比較する
    assert output_data["grna_acceptor"].tolist() == expected_output["grna_acceptor"].tolist()
    assert output_data["grna_donor"].tolist() == expected_output["grna_donor"].tolist()

def test_extract_sgrna_features():
    sgrna_list = [
        SgrnaInfo(
            target_sequence="AAA",
            actual_sequence="UUU",
            start_in_sequence=1,
            end_in_sequence=21,
            target_pos_in_sgrna=5,
            overlap_between_cds_and_editing_window=3,
            possible_unintended_edited_base_count=1,
        ),
        SgrnaInfo(
            target_sequence="CCC",
            actual_sequence="GGG",
            start_in_sequence=2,
            end_in_sequence=22,
            target_pos_in_sgrna=6,
            overlap_between_cds_and_editing_window=4,
            possible_unintended_edited_base_count=2,
        ),
    ]
    result = extract_sgrna_features(sgrna_list)
    assert result[0] == ["AAA", "CCC"]
    assert result[1] == ["UUU", "GGG"]
    assert result[2] == [1, 2]
    assert result[3] == [21, 22]
    assert result[4] == [5, 6]
    assert result[5] == [3, 4]
    assert result[6] == [1, 2]

def test_extract_acceptor_sgrna_features_empty():
    result = extract_sgrna_features([])
    assert result == ([], [], [], [], [], [], [])

def test_organize_target_exon_df_with_grna_sequence():
    acceptor_list = [
        SgrnaInfo(
            target_sequence="AAA",
            actual_sequence="UUU",
            start_in_sequence=1,
            end_in_sequence=21,
            target_pos_in_sgrna=5,
            overlap_between_cds_and_editing_window=3,
            possible_unintended_edited_base_count=1,
        )
    ]
    donor_list = [
        SgrnaInfo(
            target_sequence="CCC",
            actual_sequence="GGG",
            start_in_sequence=2,
            end_in_sequence=22,
            target_pos_in_sgrna=6,
            overlap_between_cds_and_editing_window=4,
            possible_unintended_edited_base_count=2,
        )
    ]
    df = pd.DataFrame({
        "grna_acceptor": [acceptor_list],
        "grna_donor": [donor_list],
        "other_col": [123],
    })

    result = organize_target_exon_df_with_grna_sequence(df)

    # acceptor列
    assert result["acceptor_sgrna_target_sequence"][0] == ["AAA"]
    assert result["acceptor_sgrna_actual_sequence"][0] == ["UUU"]
    assert result["acceptor_sgrna_start_in_sequence"][0] == [1]
    assert result["acceptor_sgrna_end_in_sequence"][0] == [21]
    assert result["acceptor_sgrna_target_pos_in_sgrna"][0] == [5]
    assert result["acceptor_sgrna_overlap_between_cds_and_editing_window"][0] == [3]
    assert result["acceptor_sgrna_possible_unintended_edited_base_count"][0] == [1]
    # donor列
    assert result["donor_sgrna_target_sequence"][0] == ["CCC"]
    assert result["donor_sgrna_actual_sequence"][0] == ["GGG"]
    assert result["donor_sgrna_start_in_sequence"][0] == [2]
    assert result["donor_sgrna_end_in_sequence"][0] == [22]
    assert result["donor_sgrna_target_pos_in_sgrna"][0] == [6]
    assert result["donor_sgrna_overlap_between_cds_and_editing_window"][0] == [4]
    assert result["donor_sgrna_possible_unintended_edited_base_count"][0] == [2]
    # 他の列が残っていること
    assert result["other_col"][0] == 123

def test_organize_target_exon_df_with_grna_sequence_empty():
    df = pd.DataFrame({
        "grna_acceptor": [[]],
        "grna_donor": [[]],
    })
    result = organize_target_exon_df_with_grna_sequence(df)
    for col in result.columns:
        assert result[col][0] == []

def test_modify_sgrna_start_end_position_to_position_in_chromosome():
    df = pd.DataFrame({
        "acceptor_sgrna_start_in_sequence": [[1, 5]],
        "acceptor_sgrna_end_in_sequence": [[21, 25]],
        "donor_sgrna_start_in_sequence": [[2, 6]],
        "donor_sgrna_end_in_sequence": [[22, 26]],
        "chromStart_acceptor": [100],
        "chromEnd_acceptor": [150],
        "chromStart_donor": [200],
        "chromEnd_donor": [250],
        "other_col": [123]  # 他の列も含める
    })

    result = modify_sgrna_start_end_position_to_position_in_chromosome(df)

    # ゲノム座標が正しく計算されているか
    assert result["acceptor_sgrna_start_in_genome"][0] == [101, 105]
    assert result["acceptor_sgrna_end_in_genome"][0] == [121, 125]
    assert result["donor_sgrna_start_in_genome"][0] == [202, 206]
    assert result["donor_sgrna_end_in_genome"][0] == [222, 226]
    # 他の列が残っていること
    assert result["other_col"][0] == 123

def test_modify_sgrna_start_end_position_to_position_in_chromosome_empty():
    df = pd.DataFrame({
        "acceptor_sgrna_start_in_sequence": [[]],
        "acceptor_sgrna_end_in_sequence": [[]],
        "donor_sgrna_start_in_sequence": [[]],
        "donor_sgrna_end_in_sequence": [[]],
        "chromStart_acceptor": [100],
        "chromEnd_acceptor": [150],
        "chromStart_donor": [200],
        "chromEnd_donor": [250],
    })
    result = modify_sgrna_start_end_position_to_position_in_chromosome(df)
    assert result["acceptor_sgrna_start_in_genome"][0] == []
    assert result["acceptor_sgrna_end_in_genome"][0] == []
    assert result["donor_sgrna_start_in_genome"][0] == []
    assert result["donor_sgrna_end_in_genome"][0] == []
