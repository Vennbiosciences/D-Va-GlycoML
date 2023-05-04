#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
#####################################################################################################
This script tests the glcyan de novo sequencing.
By using glycan de novo sequencing, a set of linearized glycan will output as candidates.

Created on 25 February 2022 for the unit test using pytest.
#####################################################################################################
"""
__author__ = 'ZLiang'

import pytest
from spectral_library.de_novo_sequencing import DeNovoSequencing

def test_count_total():
    # Initialization
    denovor = DeNovoSequencing()

    # Input wrong types
    with pytest.raises(ValueError):
        denovor.count_total(12040100)

    # Test Case 1
    monosaccharide_number_1 = "081204010000"
    result_1 = denovor.count_total(monosaccharide_number_1)
    # Use assert to judge the result
    assert result_1 == 25

    # Test Case 2
    monosaccharide_number_2 = "121204030102"
    result_2 = denovor.count_total(monosaccharide_number_2)
    # Use assert to judge the result
    assert result_2 == 34


def test_glycan_list_to_str():
    # Initialization
    denovor = DeNovoSequencing()

    # Input wrong types
    with pytest.raises(ValueError):
        denovor.glycan_list_to_str(12040100)

    # Test Case 1
    glycan_list_1 = "HHNFNHHNFHHFHHFHHNFNHHNFNHHNFN"
    result_1 = denovor.glycan_list_to_str(glycan_list_1)
    # Use assert to check the result
    assert result_1 == "1409070000"

    # Test Case 2
    glycan_list_2 = "HHHHHHHHHHNNNNNNNNAA"
    result_2 = denovor.glycan_list_to_str(glycan_list_2)
    # Use assert to check the result
    assert result_2 == "1008000200"

    # Test Case 3
    glycan_list_3 = "HHFHAHAAHHGGNNNFNNFNFNNFNFNF"
    result_3 = denovor.glycan_list_to_str(glycan_list_3)
    # Use assert to check the result
    assert result_3 == "0610070302"


def test_glycan_contain():
    # Initialization
    denovor = DeNovoSequencing()

    # Input wrong types
    with pytest.raises(ValueError):
        denovor.glycan_contain(12040100, "0902010000")
    with pytest.raises(ValueError):
        denovor.glycan_contain("0902010000", 12040100)

    # Test Case 1
    glycan_former_1 =  "0902010000"
    glycan_current_1 = "1003010000"
    result_1 = denovor.glycan_contain(glycan_former_1, glycan_current_1)
    # Use assert to check the result
    assert result_1 == True

    # Test Case 2
    glycan_former_2 =  "0912010000"
    glycan_current_2 = "1509010000"
    result_2 = denovor.glycan_contain(glycan_former_2, glycan_current_2)
    # Use assert to check the result
    assert result_2 == False

    # Test Case 3
    glycan_former_3 =  "0000000000"
    glycan_current_3 = "1529111111"
    result_3 = denovor.glycan_contain(glycan_former_3, glycan_current_3)
    # Use assert to check the result
    assert result_3 == True

    # Test Case 4
    glycan_former_4 =  "1529111111"
    glycan_current_4 = "1529111111"
    result_4 = denovor.glycan_contain(glycan_former_4, glycan_current_4)
    # Use assert to check the result
    assert result_4 == True

    # Test Case 5
    glycan_former_5 =  "0000000000"
    glycan_current_5 = "0000000000"
    result_5 = denovor.glycan_contain(glycan_former_5, glycan_current_5)
    # Use assert to check the result
    assert result_5 == True

    # Test Case 6
    glycan_former_6 =  "0000000010"
    glycan_current_6 = "0101010100"
    result_6 = denovor.glycan_contain(glycan_former_6, glycan_current_6)
    # Use assert to check the result
    assert result_6 == False


def test_glycan_minus():
    # Initialization
    denovor = DeNovoSequencing()

    # Input wrong types
    with pytest.raises(ValueError):
        denovor.glycan_minus(12040100, "0902010000")
    with pytest.raises(ValueError):
        denovor.glycan_minus("0902010000", 12040100)

    # Test Case 1
    glycan_former_1 =  "0902010000"
    glycan_current_1 = "1003010000"
    result_1 = denovor.glycan_minus(glycan_former_1, glycan_current_1)
    # Use assert to check the result
    assert result_1 == "0101000000"

    # Test Case 2
    glycan_former_2 =  "0912010000"
    glycan_current_2 = "1513011111"
    result_2 = denovor.glycan_minus(glycan_former_2, glycan_current_2)
    # Use assert to check the result
    assert result_2 == "0601001111"

    # Test Case 3
    glycan_former_3 =  "0000000000"
    glycan_current_3 = "1529111111"
    result_3 = denovor.glycan_minus(glycan_former_3, glycan_current_3)
    # Use assert to check the result
    assert result_3 == "1529111111"

    # Test Case 4
    glycan_former_4 =  "1529111111"
    glycan_current_4 = "1529111111"
    result_4 = denovor.glycan_minus(glycan_former_4, glycan_current_4)
    # Use assert to check the result
    assert result_4 == "0000000000"

    # Test Case 5
    glycan_former_5 =  "0000000000"
    glycan_current_5 = "0000000000"
    result_5 = denovor.glycan_minus(glycan_former_5, glycan_current_5)
    # Use assert to check the result
    assert result_5 == "0000000000"

    # Test Case 6
    glycan_former_6 =  "0101010100"
    glycan_current_6 = "2121212120"
    result_6 = denovor.glycan_minus(glycan_former_6, glycan_current_6)
    # Use assert to check the result
    assert result_6 == "2020202020"


def test_glycan_str_to_list():
    # Initialization
    denovor = DeNovoSequencing()

    # Input wrong types
    with pytest.raises(ValueError):
        denovor.glycan_str_to_list(12040100)

    # Test Case 1
    glycan_str_1 = "0912010000"
    result_1 = denovor.glycan_str_to_list(glycan_str_1)
    # Use assert to check the result
    assert result_1 == "HHHHHHHHHNNNNNNNNNNNNF"

    # Test Case 2
    glycan_str_2 = "1509010001"
    result_2 = denovor.glycan_str_to_list(glycan_str_2)
    # Use assert to check the result
    assert result_2 == "HHHHHHHHHHHHHHHNNNNNNNNNFG"

    # Test Case 3
    glycan_str_3 = "0000000000"
    result_3 = denovor.glycan_str_to_list(glycan_str_3)
    # Use assert to check the result
    assert result_3 == ""

    # Test Case 4
    glycan_str_4 =  "1529111111"
    result_4 = denovor.glycan_str_to_list(glycan_str_4)
    # Use assert to check the result
    assert result_4 == "HHHHHHHHHHHHHHHNNNNNNNNNNNNNNNNNNNNNNNNNNNNNFFFFFFFFFFFAAAAAAAAAAAGGGGGGGGGGG"

    # Test Case 5
    glycan_str_5 =  "0101010101"
    result_5 = denovor.glycan_str_to_list(glycan_str_5)
    # Use assert to check the result
    assert result_5 == "HNFAG"

    # Test Case 6
    glycan_str_6 =  "1010101010"
    result_6 = denovor.glycan_str_to_list(glycan_str_6)
    # Use assert to check the result
    assert result_6 == "HHHHHHHHHHNNNNNNNNNNFFFFFFFFFFAAAAAAAAAAGGGGGGGGGG"

    # Test Case 7
    glycan_str_7 =  "1111010101"
    result_7 = denovor.glycan_str_to_list(glycan_str_7)
    # Use assert to check the result
    assert result_7 == "HHHHHHHHHHHNNNNNNNNNNNFAG"


def test_de_novo():
    # Initialization
    denovor = DeNovoSequencing()

    # Input wrong types
    with pytest.raises(ValueError):
        denovor.de_novo(12040100, {'0001000000': 178721.8, '0002000000': 24745.4})
    with pytest.raises(ValueError):
        denovor.de_novo("0504000002", "0504000002")

    # Test Case 1
    peptide_1 = "JSSTQFEVK"
    glycan_1 = "HHHHHHHHHHNN"
    glycan_precursor_1 = "1002000000"
    glycan_intensity_dic_1 = {'0002000000': 4459.8, '0102000000': 14790.1, '0202000000': 19117.6, '0302000000': 2415.0,
                              '0402000000': 21476.3, '0001000000': 85926.2, '0602000000': 4881.8, '0702000000': 2660.0,
                              '0502000000': 12768.6, '0802000000': 4004.8}
    result_1 = denovor.de_novo(glycan_precursor_1, glycan_intensity_dic_1)
    # Use assert to check the result
    assert result_1 == [('NNHHHHHHHHHH', 172500.2)]

    # Test Case 2
    peptide_2 = "VVLHPJYSQVDIGLIK"
    glycan_2 = "HHHHHHHHHHNNNNNNNNFF"
    glycan_precursor_2 = "1008020000"
    glycan_intensity_dic_2 = {'0305010000': 7025.1, '0102000000': 221255.5, '0203000000': 9707.5, '0001000000': 164728.0,
                              '0303010000': 1166.9, '0403010000': 103254.2, '0002000000': 28635.7, '0707010000': 3439.2,
                              '0807010000': 5226.1, '0202000000': 661620.3, '0302000000': 2493.9, '0805000000': 3161.1,
                              '0606010000': 2111.2, '0507010000': 2121.4, '0203010000': 1748.1, '0303000000': 122095.3}
    result_2 = denovor.de_novo(glycan_precursor_2, glycan_intensity_dic_2)
    # Use assert to check the result
    assert result_2 == [('NNHHNHFNNHHNNHHHHHNF', 1227021.0), ('NNHHNHFNNHHHNHNHHHNF', 1227010.8),
                        ('NNHHHNFNNHHNNHHHHHNF', 1219807.4), ('NNHHHNFNNHHHNHNHHHNF', 1219797.2),
                        ('NNHHNFHNNHHNNHHHHHNF', 1106673.8), ('NNHHNFHNNHHHNHNHHHNF', 1106663.6)]

    # Test Case 3
    peptide_3 = "GTAGNALMDGASQLTGEJR"
    glycan_3 = "HHHHHHNNNNNNNNNN"
    glycan_precursor_3 = "0610000000"
    glycan_intensity_dic_3 = {'0405000000': 3821.8, '0001000000': 849741.3, '0002000000': 323886.8,
                              '0102000000': 1878235.8, '0203000000': 2047.0}
    result_3 = denovor.de_novo(glycan_precursor_3, glycan_intensity_dic_3)
    # Use assert to check the result
    assert result_3 == [('NNHHNHHNNHHNNNNN', 3057732.7)]

    # Test Case 4
    peptide_4 = "NLFLJHSEJATAK"
    glycan_4 = "HHHHHHHHHHNNNNNNNNAA"
    glycan_precursor_4 =  "1008000200"
    glycan_intensity_dic_4 = {'0001000000': 16685.6, '0002000000': 26543.0, '0202000000': 1137.5, '0103000000': 6672.6,
                              '0203000000': 3837.2, '0204000000': 1080.6, '0304000000': 4370.1, '0404000000': 3547.2,
                              '0907000100': 21704.1, '0405000000': 2264.8, '1008000000': 3106.1, '0908000100': 1285.3,
                              '0505000000': 8320.1, '1008000100': 6415.6, '0605000000': 1467.8, '0506000000': 1652.4,
                              '1008000200': 5649.7, '0407000000': 13622.1, '0505000100': 11062.8, '0706000000': 5844.5,
                              '0806000000': 3359.1, '0707000000': 1889.9, '0807000000': 8296.6, '0806000100': 3876.6,
                              '0907000000': 19080.1, '0807000100': 7476.3, '0907000200': 2312.4}
    result_4 = denovor.de_novo(glycan_precursor_4, glycan_intensity_dic_4)
    # Use assert to check the result
    assert result_4 == [('NNHNHNHHNHNHHHNHANHA', 146608.6), ('NNHNHNHHNHHHNHNHANHA', 146424.0),
                        ('NNHNHNHHNHNHHNHHANHA', 145139.4), ('NNHNHNHHNHHHNNHHANHA', 144954.8),
                        ('NNHHNNHHNHNHHHNHANHA', 141073.5), ('NNHHNNHHNHHHNHNHANHA', 140888.9),
                        ('NNHHNNHHNHNHHNHHANHA', 139604.3), ('NNHHNNHHNHHHNNHHANHA', 139419.7),
                        ('NNHNHNHHNHNHHHNAHNHA', 135004.8), ('NNHNHNHHNHHHNHNAHNHA', 134820.2),
                        ('NNHNHNHHNHNHHNHAHNHA', 133535.6), ('NNHNHNHHNHHHNNHAHNHA', 133351.0),
                        ('NNHNHNHHNHNHHHANHNHA', 130584.8), ('NNHNHNHHNHHHNHANHNHA', 130400.2),
                        ('NNHHNNHHNHNHHHNAHNHA', 129469.7), ('NNHHNNHHNHHHNHNAHNHA', 129285.1),
                        ('NNHHNNHHNHNHHNHAHNHA', 128000.5), ('NNHHNNHHNHHHNNHAHNHA', 127815.9),
                        ('NNHNHNHHNHNHHHNHHNAA', 126725.3), ('NNHNHNHHNHHHNHNHHNAA', 126540.7),
                        ('NNHNHNHHNHNHHNHHHNAA', 125256.1), ('NNHNHNHHNHHHNNHHHNAA', 125071.5),
                        ('NNHHNNHHNHNHHHANHNHA', 125049.7), ('NNHHNNHHNHHHNHANHNHA', 124865.1),
                        ('NNHHNNHHNHNHHHNHHNAA', 121190.2), ('NNHHNNHHNHHHNHNHHNAA', 121005.6),
                        ('NNHHNNHHNHNHHNHHHNAA', 119721.0), ('NNHHNNHHNHHHNNHHHNAA', 119536.4)]

    # Test Case 5
    peptide_5 = "NLFLJHSEJATAK"
    glycan_5 = "HHHHHHHHHHNNNNNNNNA"
    glycan_precursor_5 = "1008000100"
    glycan_intensity_dic_5 = {'0001000000': 16193.6, '0002000000': 17740.9, '0102000000': 1247.5, '0203000000': 3214.0,
                              '0404000000': 5335.0, '0305000000': 1276.9, '0908000000': 3783.4, '0907000100': 21201.0,
                              '0405000000': 2343.0, '1008000000': 6672.5, '0908000100': 2974.4, '0505000000': 5465.3,
                              '1008000100': 21644.1, '0605000000': 2356.8, '0506000000': 2210.6, '0505000100': 2846.3,
                              '0705000000': 1561.4, '0606000000': 1378.0, '0506000100': 978.0, '0706000000': 14191.9,
                              '0607000000': 1505.8, '0606000100': 2450.8, '0806000000': 7586.2, '0707000000': 2917.7,
                              '0706000100': 2062.3, '0807000000': 14267.8, '0806000100': 3318.2, '0907000000': 36784.5,
                              '0807000100': 5004.5}
    result_5 = denovor.de_novo(glycan_precursor_5, glycan_intensity_dic_5)
    # Use assert to check the result
    assert result_5 == [('NNHHNHHNNHHHNHNHANH', 174107.4), ('NNHHNHHNNHHNHHNHANH', 173924.0),
                        ('NNHHNHHNNHNHHHNHANH', 173777.8), ('NNHHNHNNHHHHNHNHANH', 170049.3),
                        ('NNHHNHNNHHHNHHNHANH', 169865.9), ('NNHHNHNNHHNHHHNHANH', 169719.7),
                        ('NNHHNHHNNHHHNNHHANH', 169438.9), ('NNHHNHHNNHHNHNHHANH', 169255.5),
                        ('NNHHNHHNNHNHHNHHANH', 169109.3), ('NNHHNHNNHHHHNNHHANH', 165380.8),
                        ('NNHHNHNNHHHNHNHHANH', 165197.4), ('NNHHNHNNHHNHHNHHANH', 165051.2),
                        ('NNHHNHHNNHHHNHNHNHA', 160387.9), ('NNHHNHHNNHHNHHNHNHA', 160204.5),
                        ('NNHHNHHNNHNHHHNHNHA', 160058.3), ('NNHHNHHNNHHHNHNHNAH', 156689.8),
                        ('NNHHNHHNNHHNNHHHANH', 156569.4), ('NNHHNHHNNHHNHHNHNAH', 156506.4),
                        ('NNHHNHHNNHNHNHHHANH', 156423.2), ('NNHHNHHNNHNHHHNHNAH', 156360.2),
                        ('NNHHNHNNHHHHNHNHNHA', 156329.8), ('NNHHNHNNHHHNHHNHNHA', 156146.4),
                        ('NNHHNHNNHHNHHHNHNHA', 156000.2), ('NNHHNHHNNHHHNNHHNHA', 155719.4),
                        ('NNHHNHHNNHHNHNHHNHA', 155536.0), ('NNHHNHHNNHNHHNHHNHA', 155389.8),
                        ('NNHHNHNNHHHHNHNHNAH', 152631.7), ('NNHHNHNNHHHNNHHHANH', 152511.3),
                        ('NNHHNHNNHHHNHHNHNAH', 152448.3), ('NNHHNHNNHHNHNHHHANH', 152365.1),
                        ('NNHHNHNNHHNHHHNHNAH', 152302.1), ('NNHHNHHNNHHHNNHHNAH', 152021.3),
                        ('NNHHNHHNNHHNHNHHNAH', 151837.9), ('NNHHNHHNNHNHHNHHNAH', 151691.7),
                        ('NNHHNHNNHHHHNNHHNHA', 151661.3), ('NNHHNHNNHHHNHNHHNHA', 151477.9),
                        ('NNHHNHNNHHNHHNHHNHA', 151331.7), ('NNHHNHNNHHHHNNHHNAH', 147963.2),
                        ('NNHHNHNNHHHNHNHHNAH', 147779.8), ('NNHHNHNNHHNHHNHHNAH', 147633.6),
                        ('NNHHNHHNNHHNNHHHNHA', 142849.9), ('NNHHNHHNNHNHNHHHNHA', 142703.7),
                        ('NNHHNHHNNHHHNHNAHNH', 142327.4), ('NNHHNHHNNHHNHHNAHNH', 142144.0),
                        ('NNHHNHHNNHNHHHNAHNH', 141997.8), ('NNHHNHHNNHHNNHHHNAH', 139151.8),
                        ('NNHHNHHNNHNHNHHHNAH', 139005.6), ('NNHHNHNNHHHNNHHHNHA', 138791.8),
                        ('NNHHNHNNHHNHNHHHNHA', 138645.6), ('NNHHNHNNHHHHNHNAHNH', 138269.3),
                        ('NNHHNHNNHHHNHHNAHNH', 138085.9), ('NNHHNHNNHHNHHHNAHNH', 137939.7),
                        ('NNHHNHHNNHHHNNHAHNH', 137658.9), ('NNHHNHHNNHHNHNHAHNH', 137475.5),
                        ('NNHHNHHNNHNHHNHAHNH', 137329.3), ('NNHHNHNNHHHNNHHHNAH', 135093.7),
                        ('NNHHNHNNHHNHNHHHNAH', 134947.5), ('NNHHNHNNHHHHNNHAHNH', 133600.8),
                        ('NNHHNHNNHHHNHNHAHNH', 133417.4), ('NNHHNHNNHHNHHNHAHNH', 133271.2),
                        ('NNHHNHHNNHHHNHANHNH', 131377.8), ('NNHHNHHNNHHNHHANHNH', 131194.4),
                        ('NNHHNHHNNHNHHHANHNH', 131048.2), ('NNHHNHNNHHHHNHANHNH', 127319.7),
                        ('NNHHNHNNHHHNHHANHNH', 127136.3), ('NNHHNHNNHHNHHHANHNH', 126990.1),
                        ('NNHHNHHNNHHHNAHNHNH', 125853.9), ('NNHHNHHNNHHNHAHNHNH', 125670.5),
                        ('NNHHNHHNNHNHHAHNHNH', 125524.3), ('NNHHNHHNNHHNNHHAHNH', 124789.4),
                        ('NNHHNHHNNHNHNHHAHNH', 124643.2), ('NNHHNHNNHHHHNAHNHNH', 121795.8),
                        ('NNHHNHNNHHHNHAHNHNH', 121612.4), ('NNHHNHNNHHNHHAHNHNH', 121466.2),
                        ('NNHHNHNNHHHNNHHAHNH', 120731.3), ('NNHHNHNNHHNHNHHAHNH', 120585.1),
                        ('NNHHNHHNNHANHHHNHNH', 114018.9), ('NNHHNHHNNHHNAHHNHNH', 113929.4),
                        ('NNHHNHHNNHNHAHHNHNH', 113783.2), ('NNHHNHHNNHNAHHHNHNH', 113383.2),
                        ('NNHHNHNNHHANHHHNHNH', 109960.8), ('NNHHNHNNHHHNAHHNHNH', 109871.3),
                        ('NNHHNHNNHHNHAHHNHNH', 109725.1), ('NNHHNHNNHHNAHHHNHNH', 109325.1)]

    # Test Case 6
    peptide_6 = "VVLHPJYSQVDIGLIK"
    glycan_6 = "HHHHHHHHHHHHNNNNNNNA"
    glycan_precursor_6 = "1207000100"
    glycan_intensity_dic_6 = {'0307000000': 1836.7, '0504000100': 5121.4, '0705000000': 4606.0, '0001000000': 170767.1,
                              '0002000000': 43005.5, '0907000000': 4521.5, '0102000000': 38524.9, '0507000000': 2130.6,
                              '0203000000': 1927787.0, '0303000000': 5621.4, '0403000000': 6529.3, '0303000100': 3586.3,
                              '1106000000': 4946.5, '0403000100': 9046.1}
    result_6 = denovor.de_novo(glycan_precursor_6, glycan_intensity_dic_6)
    # Use assert to check the result
    assert result_6 == [('NNHHNHHHHHNNHHHHNHNA', 2201787.7)]

    # Test Case 7
    peptide_7 = "LSLHRPALEDLLLGSEAJLTCTLTGLR"
    glycan_7 = "HHHHHNNNNA"
    glycan_precursor_7 = "0504000100"
    glycan_intensity_dic_7 = {'0001000000': 33081.3, '0504000000': 3832.5, '0202000000': 8115.3, '0302000000': 2097.4,
                              '0303000000': 12044.7, '0403000000': 32379.6, '0102000000': 2790.4}
    result_7 = denovor.de_novo(glycan_precursor_7, glycan_intensity_dic_7)
    # Use assert to check the result
    assert result_7 == [('NHNHHNHHNA', 94341.2)]

