import pandas as pd
import re
from altex_aid.sgrna_designer import (
    SgrnaInfo,
    BaseEditor,
    convert_dna_to_reversed_complement_rna,
    convert_dna_to_rna,
    reverse_complement_pam_as_regex,
    convert_pam_as_regex,
    calculate_overlap_and_unintended_edits_to_cds,
    design_sgrna,
    decide_target_base_pos_in_sequence,
    is_valid_exon_position,
    design_sgrna_for_target_exon_df,
    extract_sgrna_features,
    organize_target_exon_df_with_grna_sequence,
    convert_sgrna_start_end_position_to_position_in_chromosome,
    design_sgrna_for_base_editors
)   

pd.set_option('display.max_columns', None)  # 全てのカラムを表示するための設定

def test_convert_dna_to_reversed_complement_rna():
    input_sequence = "ATGCATGC"
    expected_output = "GCAUGCAU" # RNAにした逆相補鎖
    output_sequence = convert_dna_to_reversed_complement_rna(input_sequence)
    assert output_sequence == expected_output

def test_convert_dna_to_rna():
    input_sequence = "ATGCATGC"
    expected_output = "AUGCAUGC"  # RNAに変換した結果
    output_sequence = convert_dna_to_rna(input_sequence)
    assert output_sequence == expected_output

def test_reverse_complement_pam_as_regex():
    input_sequence = "NGG"
    expected_output = re.compile("(?=([Cc][Cc][ATGCatgc]))")
    output_sequence = reverse_complement_pam_as_regex(input_sequence)
    assert output_sequence == expected_output

def test_convert_pam_as_regex():
    input_sequence = "NGG"
    expected_output = re.compile("(?=([ATGCatgc][Gg][Gg]))")  # PAM配列を正規表現に変換した結果
    output_sequence = convert_pam_as_regex(input_sequence)
    assert output_sequence == expected_output

def test_calculate_overlap_and_unintended_edits_to_cds_cbe():
    editing_sequence = "NNNCCCCCNNNNNNNNNNNNNNNAGGGNNNNNNNNNNNNNNNNNNNNNNN"
    window_start_in_seq = 24
    window_end_in_seq = 26 # このendはinclusiveなので、24,25,26の3塩基が編集される
    cds_boundary = 25
    site_type = "acceptor"
    base_editor_type = "cbe"  # base editorのタイプを指定

    expected_output = (2,2)
    overlap, unintended_edits = calculate_overlap_and_unintended_edits_to_cds(
        editing_sequence=editing_sequence,
        window_start_in_seq=window_start_in_seq,
        window_end_in_seq=window_end_in_seq,
        cds_boundary=cds_boundary,
        site_type=site_type,
        base_editor_type=base_editor_type
    )
    assert (overlap, unintended_edits) == expected_output

def test_calculate_overlap_and_unintended_edits_to_cds_abe():
    editing_sequence = "NNNCCCCCNNNNNNNNNNNNNNNAGAAANNNNNNNNNNNNNNNNNNNNNNNN"
    window_start_in_seq = 23
    window_end_in_seq = 25
    cds_boundary = 25
    site_type = "acceptor"
    base_editor_type = "abe"  # base editorのタイプを指定

    expected_output = (1,1)
    overlap, unintended_edits = calculate_overlap_and_unintended_edits_to_cds(
        editing_sequence=editing_sequence,
        window_start_in_seq=window_start_in_seq,
        window_end_in_seq=window_end_in_seq,
        cds_boundary=cds_boundary,
        site_type=site_type,
        base_editor_type=base_editor_type
    )
    assert (overlap, unintended_edits) == expected_output

def test_design_sgrna_cbe_acceptor():
    editing_sequence = "NNNCCCCCNNNNNNNNNNNNNNNAGGGNNNNNNNNNNNNNNNNNNNNNNN"
    # これは+ strandのエキソンで、SAの配列。
    pam_sequence = "NGG"  # 例としてNGGを使用
    pam_regex = convert_pam_as_regex(pam_sequence)
    reversed_pam_regex = reverse_complement_pam_as_regex(pam_sequence)  # PAMはNGGなので、逆相補化してCCNとする
    editing_window_start_in_grna = 17
    editing_window_end_in_grna = 19
    target_base_pos_in_sequence = 24 # acceptorなら24番目のG, donorなら25番目のGが編集ターゲット
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
        pam_regex=pam_regex,
        reversed_pam_regex=reversed_pam_regex,
        editing_window_start_in_grna=editing_window_start_in_grna,
        editing_window_end_in_grna=editing_window_end_in_grna,
        target_base_pos_in_sequence=target_base_pos_in_sequence,
        cds_boundary=cds_boundary,
        base_editor_type="cbe",
        site_type=site_type
    )
    print(output_data)
    assert output_data == expected_output

def test_design_sgrna_cbe_donor():
    editing_sequence = "NNNNCCCCCNNNNNNNNNNNNNNGGGTNNNNNNNNNNNNNNNNNNNNN"
    # これは+ strandのエキソンで、SDの配列。
    pam_sequence = "NGG"  # 例としてNGGを使用
    pam_regex = convert_pam_as_regex(pam_sequence)
    reversed_pam_regex = reverse_complement_pam_as_regex(pam_sequence)  # PAMはNGGなので、逆相補化してCCNとする
    editing_window_start_in_grna = 17
    editing_window_end_in_grna = 19
    target_base_pos_in_sequence = 25 # acceptorなら24番目のG, donorなら25番目のGが編集ターゲット
    cds_boundary = 24
    site_type = "donor"

    expected_output = [
        SgrnaInfo(
            target_sequence="CCNNNNNNNNNNNNNNGGGT",
            actual_sequence="ACCCNNNNNNNNNNNNNNGG",
            start_in_sequence = 7,
            end_in_sequence=27,
            target_pos_in_sgrna=19,
            overlap_between_cds_and_editing_window=2,
            possible_unintended_edited_base_count=2
        ),
        SgrnaInfo(
            target_sequence="CNNNNNNNNNNNNNNGGGTN",
            actual_sequence="NACCCNNNNNNNNNNNNNNG",
            start_in_sequence = 8,
            end_in_sequence=28,
            target_pos_in_sgrna=18,
            overlap_between_cds_and_editing_window=1,
            possible_unintended_edited_base_count= 1
        ),
        SgrnaInfo(
            target_sequence="NNNNNNNNNNNNNNGGGTNN",
            actual_sequence="NNACCCNNNNNNNNNNNNNN",
            start_in_sequence = 9,
            end_in_sequence=29,
            target_pos_in_sgrna=17,
            overlap_between_cds_and_editing_window=0,
            possible_unintended_edited_base_count= 0
        ),
]
    output_data = design_sgrna(
        editing_sequence=editing_sequence,
        pam_regex=pam_regex,
        reversed_pam_regex=reversed_pam_regex,
        editing_window_start_in_grna=editing_window_start_in_grna,
        editing_window_end_in_grna=editing_window_end_in_grna,
        target_base_pos_in_sequence=target_base_pos_in_sequence,
        cds_boundary=cds_boundary,
        base_editor_type="cbe",
        site_type=site_type
    )
    print(output_data)
    assert output_data == expected_output

def test_design_sgrna_abe_acceptor():
    # acceptorの場合のテスト
    editing_sequence = "NNNNNNNNNNNNNNNNNNNNNNNAGAANNNNNNNNNNNNNGGGGGNNNNN"
    pam_sequence = "NGG"
    reversed_pam_regex = reverse_complement_pam_as_regex(pam_sequence)
    pam_regex = convert_pam_as_regex(pam_sequence)
    editing_window_start_in_grna = 17
    editing_window_end_in_grna = 19
    target_base_pos_in_sequence = 23
    cds_boundary = 25
    site_type = "acceptor"


    expected_output = [
        SgrnaInfo(
            target_sequence="NNNAGAANNNNNNNNNNNNN",
            actual_sequence="NNNAGAANNNNNNNNNNNNN",
            start_in_sequence= 20,
            end_in_sequence= 40,
            target_pos_in_sgrna=17,
            overlap_between_cds_and_editing_window=0,
            possible_unintended_edited_base_count=0
        ),
        SgrnaInfo(
            target_sequence="NNAGAANNNNNNNNNNNNNG",
            actual_sequence="NNAGAANNNNNNNNNNNNNG",
            start_in_sequence= 21,
            end_in_sequence= 41,
            target_pos_in_sgrna=18,
            overlap_between_cds_and_editing_window=0,
            possible_unintended_edited_base_count=0
        ),
        SgrnaInfo(
            target_sequence="NAGAANNNNNNNNNNNNNGG",
            actual_sequence="NAGAANNNNNNNNNNNNNGG",
            start_in_sequence= 22,
            end_in_sequence= 42,
            target_pos_in_sgrna=19,
            overlap_between_cds_and_editing_window=1,
            possible_unintended_edited_base_count=1
        ),

    ]
    output_data = design_sgrna(
        editing_sequence=editing_sequence,
        reversed_pam_regex=reversed_pam_regex,
        pam_regex=pam_regex,
        editing_window_start_in_grna=editing_window_start_in_grna,
        editing_window_end_in_grna=editing_window_end_in_grna,
        target_base_pos_in_sequence=target_base_pos_in_sequence,
        cds_boundary=cds_boundary,
        base_editor_type="abe",
        site_type=site_type
    )
    print(output_data)
    assert output_data == expected_output

def test_design_sgrna_abe_donor():
    editing_sequence = "NNNNNCCCCCNNNNNNNNNNNNNTTGTNNNNNNNNNNNNNNNNNNNNN"
    # これは+ strandのエキソンで、SDの配列。
    pam_sequence = "NGG"  # 例としてNGGを使用
    reversed_pam_regex = reverse_complement_pam_as_regex(pam_sequence)  # PAMはNGGなので、逆相補化してCCNとする
    pam_regex = convert_pam_as_regex(pam_sequence)
    editing_window_start_in_grna = 17
    editing_window_end_in_grna = 19
    target_base_pos_in_sequence = 26 # acceptorなら24番目のG, donorなら25番目のGが編集ターゲット
    cds_boundary = 24
    site_type = "donor"

    expected_output = [
        SgrnaInfo(
            target_sequence="CCNNNNNNNNNNNNNTTGTN",
            actual_sequence="NACAANNNNNNNNNNNNNGG",
            start_in_sequence = 8,
            end_in_sequence=28,
            target_pos_in_sgrna=19,
            overlap_between_cds_and_editing_window=1,
            possible_unintended_edited_base_count=1
        ),
        SgrnaInfo(
            target_sequence="CNNNNNNNNNNNNNTTGTNN",
            actual_sequence="NNACAANNNNNNNNNNNNNG",
            start_in_sequence = 9,
            end_in_sequence=29,
            target_pos_in_sgrna=18,
            overlap_between_cds_and_editing_window=0,
            possible_unintended_edited_base_count= 0
        ),
        SgrnaInfo(
            target_sequence="NNNNNNNNNNNNNTTGTNNN",
            actual_sequence="NNNACAANNNNNNNNNNNNN",
            start_in_sequence = 10,
            end_in_sequence=30,
            target_pos_in_sgrna=17,
            overlap_between_cds_and_editing_window=0,
            possible_unintended_edited_base_count= 0
        ),
]
    output_data = design_sgrna(
        editing_sequence=editing_sequence,
        reversed_pam_regex=reversed_pam_regex,
        pam_regex=pam_regex,
        editing_window_start_in_grna=editing_window_start_in_grna,
        editing_window_end_in_grna=editing_window_end_in_grna,
        target_base_pos_in_sequence=target_base_pos_in_sequence,
        base_editor_type="abe",
        cds_boundary=cds_boundary,
        site_type=site_type
    )
    print(output_data)
    assert output_data == expected_output

def test_decide_target_base_pos_in_sequence():
    # テストケース: CBE + acceptor
    base_editor_type = "cbe"
    site_type = "acceptor"
    expected_output = 24
    assert decide_target_base_pos_in_sequence(base_editor_type, site_type) == expected_output

    # テストケース: CBE + donor
    base_editor_type = "cbe"
    site_type = "donor"
    expected_output = 25
    assert decide_target_base_pos_in_sequence(base_editor_type, site_type) == expected_output

    # テストケース: ABE + acceptor
    base_editor_type = "abe"
    site_type = "acceptor"
    expected_output = 23
    assert decide_target_base_pos_in_sequence(base_editor_type, site_type) == expected_output

    # テストケース: ABE + donor
    base_editor_type = "abe"
    site_type = "donor"
    expected_output = 26
    assert decide_target_base_pos_in_sequence(base_editor_type, site_type) == expected_output


def test_is_valid_exon_position_acceptor():
    # acceptor の場合のテスト
    assert is_valid_exon_position("internal", "acceptor") is True
    assert is_valid_exon_position("last", "acceptor") is True
    assert is_valid_exon_position("first", "acceptor") is False

def test_is_valid_exon_position_donor():
    # donor の場合のテスト
    assert is_valid_exon_position("internal", "donor") is True
    assert is_valid_exon_position("first", "donor") is True
    assert is_valid_exon_position("last", "donor") is False

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
        "acceptor_exon_intron_boundary_±25bp_sequence":["NNNNCCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNN","NNNNCCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNN"],
        "chromStart_donor":[175, 375],
        "chromEnd_donor":[225, 425],
        "donor_exon_intron_boundary_±25bp_sequence":["NNNNNCCCNNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNNNNNNNNNNNN","NNNNNCCCNNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNNNNNNNNNNNN"]
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
        "acceptor_exon_intron_boundary_±25bp_sequence":["NNNNCCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNN","NNNNCCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNN"],
        "chromStart_donor": [175, 375],
        "chromEnd_donor": [225, 425],
        "donor_exon_intron_boundary_±25bp_sequence":["NNNNNCCCNNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNNNNNNNNNNNN","NNNNNCCCNNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNNNNNNNNNNNN"],
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
        base_editor_type="cbe",
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

def test_convert_sgrna_start_end_position_to_position_in_chromosome():
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

    result = convert_sgrna_start_end_position_to_position_in_chromosome(df)

    # ゲノム座標が正しく計算されているか
    assert result["acceptor_sgrna_start_in_genome"][0] == [101, 105]
    assert result["acceptor_sgrna_end_in_genome"][0] == [121, 125]
    assert result["donor_sgrna_start_in_genome"][0] == [202, 206]
    assert result["donor_sgrna_end_in_genome"][0] == [222, 226]
    # 他の列が残っていること
    assert result["other_col"][0] == 123

def test_convert_sgrna_start_end_position_to_position_in_chromosome_empty():
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
    result = convert_sgrna_start_end_position_to_position_in_chromosome(df)
    assert result["acceptor_sgrna_start_in_genome"][0] == []
    assert result["acceptor_sgrna_end_in_genome"][0] == []
    assert result["donor_sgrna_start_in_genome"][0] == []
    assert result["donor_sgrna_end_in_genome"][0] == []