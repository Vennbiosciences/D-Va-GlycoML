#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read MS/MS spectral library data from an MSP format file

Created on 17 September 2021.
Modified on 16 February 2022 for five monosaccharides.
Modified on 21 February 2022 for double-digit monosaccharides.
Modified on 24 February 2022 for combining 'I' and 'L' into 'X', to support the fifth monosaccharide.
Modified on 28 February 2022 for handling the fixed modification for amino acid 'C', such as Carbamidomethyl[C].
Modified on 04 March 2022, for the filters of large peptides and glycans.
Modified on 08 April 2022, for the handling of N and O linked glycopeptides.
Modified on 21 April 2022, for only outputting csv files with only one denovo candidate.
Modified on 26 April 2022, for outputting csv files with top N denovo candidates for CLI.
Modified on 11 August 2022, filter out those ions with charge > 4.
###################################################################################################################
"""
__author__ = 'ZLiang'


import csv
import re
from stack_queue import Stack
from de_novo_sequencing import DeNovoSequencing
import argparse
from pathlib import Path
import time


"""
The following is a sample spectrum for an msp file for five monosaccharides.

Name: 202110_Palleon_Cyno_In_vitro-20211027_Palleon_A_HILIC.3200.3200.5.0.dta
MW: 752.10299
Comment: CHARGE=5+      RTINSECONDS=351.1281    PEPTIDE=HPHJJSSDLHPHK   MODIFICATIONS=nan
GLYCAN(H,N,F,A,G)=5,4,0,0,2     GLYCANS=H(5)N(4)G(2)
Num Peaks: 29
235.119 9667.3  b2+1
284.172 834.2   y2+1
381.224 5364.4  y3+1
631.368 1173.7  y5+1
643.303 1563.7  y11+2
684.824 795.1   y$11+2
695.045 1616.6  Y-H(4)N(3)+4
744.844 955.5   y11-N(1)+2
750.665 1227.9  Y-H(2)N(2)+3
771.821 768.5   Y-H(4)N(3)G(1)+4
861.899 4465.2  Y-N(1)+2
863.096 1567.3  Y-H(5)N(4)G(1)+4
872.373 3981.9  Y-H(3)N(3)+3

"""

beginning_modifier_regex = r"^[a-z]\[\d+\]"
middle_modifier_regex = r"[A-Z]\[\d+\]"
position_regex = r"^[by](?P<position>\d+)"
charge_regex = r"\+(?P<charge>\d)"
backbone_ions = "by"
glyco_ions = "Y"
precursor_ions = "M"
intact_glycopeptide = "M"
monosaccharide_list = "HNFAG"

# Assume the maximum length of peptide is 32, and maximum length of glycan is 18.
MAX_PEPTIDE_LENGTH = 32
MAX_GLYCAN_LENGTH = 18

# Headers for the csv output
csv_output_rows = [
    'glycopeptide',
    'charge',
    'precursor_mass',
    'retention_time',
    'mzs',
    'intensities',
    'ions',
    'positions',
    'ion_charges',
]

# Using special characters to represent four monosaccharides
monosaccharide_component_replacements = {
    "H": "!",
    "N": "@",
    "F": "#",
    "A": "$",
    "G": "%",
}

# Combine 20 amino acids and 5 monosaccharides together to generate 25 codes, the total number is 26 after plus zero.
# An additional amino acid is 'J', which represents the glycosylated cite for 'N'.
# Ignore 'X' which is used to represent both 'I' and 'L'.
# Combine 'I' and 'L' into 'X', save the room for the fifth monosaccharide.
# Using 'I' and 'L' to represent monosaccharides, and 'Z' to represent empty one.
# The original one:  amino_acid_codes = "ACDEFGHIJKLMNPQRSTVWY"
amino_acid_codes = "ACDEFGHJKMNPQRSTVWXY"
monosaccharide_codes = "".join(monosaccharide_component_replacements.values())
zero_code = "Z"
amino_acid_monosaccharide_zero_codes = amino_acid_codes + monosaccharide_codes + zero_code
#print(amino_acid_monosaccharide_zero_codes)

# Transfer the sequence into one hot encodings by using the code
def one_hot_encode(values, code):
    """
    One hot encoding for the values by using the code.
    :param values: A string whose letter belongs to the code, such as "XSXHRPAXEDXXXGSEAJXTCTXTGXRZZZZZ@!@!!@!!@$ZZZZZZZZ".
    :param code: A string contains all the letters for the values, such as "ACDEFGHJKMNPQRSTVWXY!@#$%Z".
    :return: One hot encoding for the string, such as
        [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ...,
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]
    """
    result = []
    for letter in values:
        letter_encoding = [1 if code_letter == letter else 0 for code_letter in code]
        try:
            assert (sum(letter_encoding) == 1)
        except AssertionError:
            print("One hot encoding failed - unexpected letter in this name")
            print(values)
            raise
        result.append(letter_encoding)
    return result

# Transfer one hot encodings back into the sequence
def reverse_one_hot_encode(vectors, code):
    """
    Get the corresponding sequence from one hot encodings, by using the code.
    :param vectors: A list of one hot encodings, such as
    [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ...,
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]
    :param code: a string contains the code for the one hot, such as "ACDEFGHJKMNPQRSTVWXY!@#$%Z".
    :return: a string converted from the vectors, by using the code, such as
        "VVXHPJYSQVDXGXXKZZZZZZZZZZZZZZZZ@@!!@!@!!@!@!$ZZZZ".
    """
    letters = []
    for vector in vectors:
        i = vector.index(1)  # get the index of the item which is 1
        letters.append(code[i])
    return "".join(letters)

# reverse_one_hot_encode(one_hot_encode(values, code), code) should be the values.

def amino_acid_monosaccharide_name_parse(peptide_name, glycan_de_novo_sequence):
    """ Combine peptide and glycan to generate one hot encoding with length of 50
    :param peptide_name: name string with amino acids.
           glycan_de_novo_sequence: glycan de novo sequence with monosaccharide.
    :return: glycopeptide_one_hot_encoding
    """
    # Assume the maximum length of peptide is 32, pad the left locations with "Z"
    left_peptide_length = MAX_PEPTIDE_LENGTH - len(peptide_name)
    peptide_pad_name = peptide_name + 'Z' * left_peptide_length
    # Replace the monosaccharides with characters
    for monosaccharide, monosaccharide_code in monosaccharide_component_replacements.items():
        # Replace monosaccharide with our special codes
        glycan_de_novo_sequence = glycan_de_novo_sequence.replace(monosaccharide, monosaccharide_code)
    # Calculate the left length of glycan
    left_glycan_length = MAX_GLYCAN_LENGTH - len(glycan_de_novo_sequence)
    glycan_pad_name = glycan_de_novo_sequence + 'Z' * left_glycan_length
    glycopeptide_name = peptide_pad_name + glycan_pad_name
    return one_hot_encode(glycopeptide_name, amino_acid_monosaccharide_zero_codes)


def reverse_amino_acid_monosaccharide_coding(vectors):
    """
    Reverse one hot encoding and character substitutions for amino acids and monosaccharide
    :param vectors: one hot encoded vectors
    :return: original amino acid name and glycan string
    """
    return reverse_one_hot_encode(vectors)


def parse_glycan(glycan):
    """
    Parse glycan text such as H(2)N(2)F(1), into a glycan dictionary {H:2, N:2, F:1}
    :param glycan: A string with [character(digit+)]+
    :return: A dictionary with {character: digit}
    """
    # Also consider the following situation for double-digit monosaccharide: "H(50)N(40)F(10)"
    stack = Stack()
    glycan_dic = {}
    for i in glycan:
        if i == ")":
            monosaccharide_number = stack.pop()
            check_ahead = stack.pop()
            # Assume there are maximum double-digit for a monosaccharide, and one character for a type
            if check_ahead != "(":
                monosaccharide_number = check_ahead + monosaccharide_number
                # the former one should be "(", pop to ignore it
                stack.pop()
            monosaccharide_number = int(monosaccharide_number)
            monosaccharide_type = stack.pop()
            glycan_dic[monosaccharide_type] = monosaccharide_number
        else:
            stack.push(i)
    return glycan_dic


def glycan_to_str(glycan_dic):
    """
    Consider double-digit monosaccharide, use two digits to represent a monosaccharide.
    Parse glycan dictionary such as {H:2, N:10, F:1} into a string "0210010000", which represents the number for "HNFAG"
    :param glycan_dic: A dictionary with {character: digit}
    :return: A string like "0210010000"
    """
    glycan_mono = "HNFAG"
    glycan_number = ['0', '0', '0', '0', '0']
    for key, value in glycan_dic.items():
        glycan_number[glycan_mono.find(key)] = str(value)
    # should add 0 to the front of glycan number if it is only one digit
    for i in range(len(glycan_number)):
        if len(glycan_number[i]) == 1:
            glycan_number[i] = '0' + glycan_number[i]

    return "".join(glycan_number)

"""
def create_annoying_string(row_data, row_id):

    Row data is dictionary like follows:
    row_dictionary = {
                    'acetyl': int(has_beginning_charge),
                    'name': name_one_hot_encoded,
                    'charge': charge,
                    'precursor': precursor_mz,
                    'mz': m_over_z,
                    'intensity': intensities,
                    'ion': ions,
                    'position': positions,
                    'neutral_loss': neutral_losses,
                    'ion_charge': ion_charges,
                    'delta': deltas,
                }
    You can read it like this using CSV dictreader.
    :return: old string.

    # Do some conversions
    acetyl = str(row_data['acetyl']) == "True"
    name_one_hot_encoded = eval(str(row_data['name']))
    charge = int(row_data['charge'])
    precursor = float(row_data['precursor'])
    mzs = eval(str(row_data['mz']))
    intensities = eval(str(row_data['intensity']))
    ions = eval(str(row_data['ion']))
    positions = eval(str(row_data['position']))
    neutral_losses = eval(str(row_data['neutral_loss']))
    ion_charges = eval(str(row_data['ion_charge']))
    deltas = eval(str(row_data['delta']))


    # Reverse the OHE
    name = reverse_amino_acid_coding(name_one_hot_encoded, charge, acetyl)

    # Get the end string:
    ends = []
    for mz, intensity, ion, position, neutral_loss, ion_charge, delta in zip(
            mzs, intensities, ions, positions, neutral_losses, ion_charges, deltas):
        ion = backbone_ions[int(ion)]
        s = f"{mz}\t{intensity}\t{ion}{position}"
        if ion_charge != 1:
            s += f"^{ion_charge}"
        s += f"/{delta}"
        ends.append(s)

    #return fName: {name}
    #        LibID: {row_id}
    #        MW: {round(charge*precursor, 4)}
    #        PrecursorMZ: {precursor}
    #        Status: Normal
    #        FullName: {name}
    #        Comment: No Comment
    #        NumPeaks: {len(intensities)}\n" +"\n".join(ends)+"\n"
    """


# read the processed information from a msp file, and write the data into a csv file and a denovo msp file.
def msp_to_csv_denovo(top_number, msp_file, csv_file, denovo_file):
    """
    Write the information of input_csv into the files of csv and denovo msp;
    :param top_number: the top number of denovo results;
    :param msp_file: the msp file to read;
    :param csv_file: the csv file to write;
    :param denovo_file: the denovo msp file to write.
    :return: the total number of samples for a csv file.
    """
    num_samples = 0
    num_dl_spectra = 0
    num_spectra = 0
    with open(msp_file) as spec_library_file, \
            open(csv_file, "w", newline='') as csv_output_file, \
            open(denovo_file, 'w') as denovo_msp_file:
        writer_csv = csv.DictWriter(csv_output_file, fieldnames=csv_output_rows)
        writer_csv.writeheader()

        while True:
            current_chunk = []
            read = False
            for line in spec_library_file:
                read = True
                # Initialization
                denovor = DeNovoSequencing()
                # Use blank line to determine the start and end of each spectrum
                if not line.strip():
                    # End of chunk - start processing

                    # Variables to write:
                    # encoded_name
                    # charge

                    # Write the processed data into a csv file.
                    precursor_mass_line = current_chunk[1]
                    try:
                        assert (precursor_mass_line.startswith("MW: "))
                    except AssertionError:
                        print("This chunk has some extra lines, or some deleted lines!")
                        print(current_chunk)
                        break

                    num_spectra += 1
                    # Check the 'MODIFICATIONS' in Comment, there are two modifications right now.
                    # Example: MODIFICATIONS = 10, Carbamidomethyl[C];7, Oxidation[M];
                    # Fixed modification: Carbamidomethyl[C]
                    # Variable modification: Oxidation[M]
                    # We ignore the variable modification and keep the fixed modification, or ignore all
                    # spectra that have the variable modifications.
                    comment, charge, retention_time, peptide, modifications, glycan, glycans = current_chunk[2].split()
                    if modifications.find("Oxidation[M]") != -1:
                        #print("Has variable modifications, ignore it!")
                        #print(peptide)
                        #print(modifications)
                        break

                    # Parse precursor charge, retention time, peptide backbone, glycan compositions
                    # Only store the number of charge which is smaller than 10, ignore the '+'
                    charge = charge.split("=")[1][:1]
                    retention_time = retention_time.split("=")[1]
                    peptide_name = peptide.split("=")[1]

                    # Add filters to remove peptides with long sequences, whose lengths are
                    # greater than the maximum thresholds.
                    # For the lengths of peptides, 0.6% are greater than 32;
                    if len(peptide_name) > MAX_PEPTIDE_LENGTH:
                        #print("Peptide length is greater than 32, ignore it!")
                        #print(peptide)
                        break

                    # Use 'X' to represent 'I' and 'L' that with the same masses,
                    # save room for the fifth monosaccharide.
                    peptide_name = peptide_name.replace('I', 'X')
                    peptide_name = peptide_name.replace('L', 'X')
                    # Assume the maximum length of peptide is 32, pad the left locations with "Z"
                    left_peptide_length = MAX_PEPTIDE_LENGTH - len(peptide_name)
                    peptide_pad_name = peptide_name + 'Z' * left_peptide_length
                    # print(peptide_pad_name)
                    glycan_compositions = glycan.split("=")[1]
                    # Parse monosaccharides from glycan compositions: GLYCAN(H,N,F,A,G)=5,4,1,1,1
                    # Should consider double-digit monosaccharides.
                    monosaccharide_numbers = []
                    if glycan_compositions == "nan":
                        monosaccharide_numbers = ['0', '0', '0', '0', '0']
                    else:
                        monosaccharide_numbers = glycan_compositions.split(",")
                    # Transfer the number of monosaccharides to glycan compositions such as "5410"
                    monosaccharides = ''
                    glycan_length = 0
                    for i in range(len(monosaccharide_numbers)):
                        number = int(monosaccharide_numbers[i])
                        for _ in range(number):
                            monosaccharides += monosaccharide_list[i]
                        glycan_length += number

                    # Add filter to remove glycans with long sequences, whose lengths are
                    # greater than the maximum threshold.
                    # For the lengths of glycans, 0.2% are greater than 18.
                    if glycan_length > MAX_GLYCAN_LENGTH:
                        #print("Glycan length is greater than 18, ignore it!")
                        #print(peptide)
                        #print(glycan_compositions)
                        break

                    precursor_mass = float(precursor_mass_line.split()[1])

                    # Check the line of "Num Peaks", if is 0, ignore the following processing:
                    number_peaks_line = current_chunk[3]
                    try:
                        assert (number_peaks_line.split(":")[1].strip() != '0')
                    except AssertionError:
                        print("This spectrum has no peaks!")
                        continue

                    m_over_z_s = []
                    intensities = []
                    ions = []
                    ion_charges = []
                    positions = []

                    glycan_intensity = {}
                    # Read MS2 ions
                    for data in current_chunk[4:]:
                        mz, intensity, annotation = data.split()
                        # Only consider maximum charge of 4 now.
                        charge_search = re.search(charge_regex, annotation)
                        ion_charge = int(charge_search.group("charge"))
                        # ignore the ions with charge 5, 6
                        if ion_charge > 4:
                            continue
                        # For b$2+1, $ means "cross ring fragment of a HexNAc", should be an 83Da, ignore it now
                        # For Y0+2, Y means "peptide backbone charge 2", would ignore it now
                        # For y12-N(1)+2, peptide fragment plus a glycan, also ignore it now.
                        if annotation[0] in backbone_ions:
                            # b, or y ions in the peptide backbone
                            # b$2+1, $ means "cross ring fragment of a HexNAc", should be an 83Da
                            # y12-N(1)+2, peptide fragment plus a glycan
                            # For b ions: b3 is 0, b$3 is 1, b3-N(1) is 2,
                            # For y ions: y3 is 3, y$3 is 4, y3-N(1) is 5.
                            ion_indicator = backbone_ions.index(annotation[0])
                            if '$' in annotation:
                                ion_indicator_code = ion_indicator * 3 + 1
                                annotation = annotation.replace('$', '')
                            elif '-' in annotation:
                                ion_indicator_code = ion_indicator * 3 + 2
                            else:
                                ion_indicator_code = ion_indicator * 3
                            position_search = re.search(position_regex, annotation)
                            position = int(position_search.group("position"))
                            # Consider padding empty amino acids to the length of 32
                            if annotation[0] == 'y':
                                position += left_peptide_length
                        elif annotation[0] in glyco_ions:
                            # intact charged peptide without glycan
                            # Y0+2=1113.098
                            if annotation[1] == '0':
                                ion_indicator_code = 6
                                position = MAX_PEPTIDE_LENGTH
                            # intact charged peptide with cross ring fragment of a HexNAc
                            # Y$+2=1154.617
                            elif annotation[1] == '$':
                                ion_indicator_code = 7
                                position = MAX_PEPTIDE_LENGTH
                            # Y ions: peptide backbone attached a fragmented glycan
                            # Summarize ion intensities for each glycan composition
                            # For GLYCAN(H,N,F,A,G)=5,4,1,0,0 such as
                            # 1018.423	7877.0	Y-H(5)N(4)F(1)+3
                            # First parse H(5)N(4)F(1), then call Stack to transfer to a dictionary
                            # {'H': 5, 'N': 4, 'F': 1}, and reverse it to 5410 by adding intensity;
                            # after that, use a dictionary such as {5410: 7877.0} to represent it;
                            # by considering multiple charges, we add all of them together.
                            # annotation[1] == '-':
                            else:
                                glycan_ion = annotation.split("-")[1].split("+")[0]
                                glycan_ion_dic = parse_glycan(glycan_ion)
                                # print(glycan_ion_dic)
                                glycan_str = glycan_to_str(glycan_ion_dic)
                                # print(glycan_str)
                                if glycan_str in glycan_intensity:
                                    sum_intensity = round(float(glycan_intensity[glycan_str]) + float(intensity), 1)
                                else:
                                    sum_intensity = round(float(intensity), 1)
                                # if sum_intensity > float(intensity):
                                # print(intensity)
                                # print(sum_intensity)
                                glycan_intensity[glycan_str] = round(float(sum_intensity), 1)
                                #print(glycan_intensity)
                                # Position is counted as the total number of glycan ions, plus peptide maximum length 32
                                ion_indicator_code = 8
                                position = denovor.count_total(glycan_str) + MAX_PEPTIDE_LENGTH
                        #elif annotation[0] in precursor_ions:
                            # M-ions: intact glycopeptide (peptide attached the whole glycan)
                            # MW: 1030.0249
                            # Comment: CHARGE=5+	RTINSECONDS=5680.6674	PEPTIDE=SVQEIQATFFYFTPJK
                            # MODIFICATIONS=nan	GLYCAN(H,N,F,A,G)=7,6,0,3,0 	GLYCANS=H(7)N(6)A(3)
                            # 1030.023	1885.5	M+5
                            #glycan_ion_dic = parse_glycan(glycans)
                            # print(glycan_ion_dic)
                            #glycan_str = glycan_to_str(glycan_ion_dic)
                            # print("M ions:", glycan_str)
                            #if glycan_str in glycan_intensity:
                                #sum_intensity = round(float(glycan_intensity[glycan_str]) + float(intensity), 1)
                            #else:
                                #sum_intensity = round(float(intensity), 1)
                            # if sum_intensity > float(intensity):
                            # print(intensity)
                            # print(sum_intensity)
                        #    glycan_intensity[glycan_str] = round(float(sum_intensity), 1)
                            # print(glycan_intensity)
                            # Position is counted as the total number of glycan ions, plus peptide maximum length 32
                        #    ion_indicator_code = 2
                        #    position = denovor.count_total(glycan_str) + max_peptide_length
                            # print("M ions position:", position)
                        #    positions.append(position)
                        else:
                            # Ignore the position, use 32 + 18 to replace it
                            # ion_indicator_code = 3
                            # positions.append(50)
                            continue

                        positions.append(position)
                        m_over_z_s.append(float(mz))
                        intensities.append(round(float(intensity), 1))
                        ions.append(ion_indicator_code)
                        ion_charges.append(ion_charge)

                    glycan_sequence = glycans.split("=")[1]
                    glycan_dic = parse_glycan(glycan_sequence)
                    glycan_str = glycan_to_str(glycan_dic)
                    # Call de novo sequencing to generate linearized glycan string
                    glycan_deno_lists = denovor.de_novo(glycan_str, glycan_intensity)
                    """
                    # Output the glycans with double-digit monosaccharides.
                    for i in range(len(monosaccharide_numbers)):
                        # for two-digit monosaccharide
                        if int(monosaccharide_numbers[i]) > 9:
                            print(peptide)
                            print(modifications)
                            print("Glycan=" + monosaccharides)
                            print("Glycan String=" + glycan_str)
                            print(glycan_intensity)
                            print("Sorted denovo candidates by intensity:")
                            print(glycan_deno_lists)
                    # for the fixed modifications: "Carbamidomethyl[C]"
                    if modifications.find("Carbamidomethyl[C]") != -1:
                        print(peptide)
                        print(modifications)
                        print("Glycan=" + monosaccharides)
                        print("Glycan String=" + glycan_str)
                        print(glycan_intensity)
                        print("Sorted denovo candidates by intensity:")
                        print(glycan_deno_lists)
                    """
                    # Assume the maximum length of glycan is max_glycan_length, pad the left locations with "Z"
                    left_glycan_length = MAX_GLYCAN_LENGTH - len(monosaccharides)
                    # The number of spectra for deep learning could have several denovo results for the samples
                    num_dl_spectra += 1
                    # Consider top 10 de novo sequence candidates if the total number greater than 10
                    # top_ten_candidates = min(len(glycan_deno_lists), 10)
                    top_one_candidate = min(len(glycan_deno_lists), top_number)
                    for i in range(top_one_candidate):
                        # glycan_deno_lists[i] is a tuple like ('NNHHNHHANHNHAA', 954586.3)
                        glycan_de_novo_sequence = glycan_deno_lists[i][0]
                        for monosaccharide, monosaccharide_code in monosaccharide_component_replacements.items():
                            # Replace monosaccharide with our special codes
                            glycan_de_novo_sequence = glycan_de_novo_sequence.replace(monosaccharide,
                                                                                      monosaccharide_code)
                        glycan_pad_name = glycan_de_novo_sequence + 'Z' * left_glycan_length
                        glycopeptide_name = peptide_pad_name + glycan_pad_name
                        # print(glycopeptide_name)
                        glycopeptide_one_hot_encoded = one_hot_encode(glycopeptide_name,
                                                                      amino_acid_monosaccharide_zero_codes)
                        row_dictionary = {
                            'glycopeptide': glycopeptide_one_hot_encoded,
                            'charge': charge,
                            'precursor_mass': precursor_mass,
                            'retention_time': retention_time,
                            'mzs': m_over_z_s,
                            'intensities': intensities,
                            'ions': ions,
                            'positions': positions,
                            'ion_charges': ion_charges,
                        }
                        assert len(m_over_z_s) == len(intensities) == len(ions) == len(positions) == len(ion_charges)
                        writer_csv.writerow(row_dictionary)
                        num_samples += 1

                    # Write the msp file and denovo results into a denovo msp file.
                    # Here the denovo msp file do not contain variable modifications.
                    denovo_msp_file.write(current_chunk[0] + '\n')
                    denovo_msp_file.write(current_chunk[1] + '\n')
                    denovo_msp_file.write(current_chunk[2] + '\t' + 'DENOVO=' + str(glycan_deno_lists) + '\n')
                    denovo_msp_file.write(current_chunk[3] + '\n')
                    for data in current_chunk[4:]:
                        denovo_msp_file.write(data + '\n')
                    denovo_msp_file.write('\n')

                    # Should use break here, to successfully end the process of one chunk, and start a new one.
                    break
                else:
                    current_chunk.append(line.strip())
            # If there is no new data, end the reading process.
            if not read:
                break
    print(f"The total number of spectra is: {num_spectra}")
    print(f"The total number of spectra for the deep learning is: {num_dl_spectra}")
    print(f"The total number of samples is: {num_samples}")
    print(f"Constructed a csv file: {csv_file.name}" )
    print(f"Constructed a denovo msp file: {denovo_file.name}")
    return num_spectra, num_dl_spectra, num_samples


# Read the path for input folder (MSP) and the top number of de novo sequencing,
# then write csv files into CSV folder, and denovo msp files into DENOVO folder.
def parse_msp_files(input_type, top_number, input_path):
    """
    Write the CSV files to the folder of "input_path"
    :param input_type: A type for the input folder, such as "N" or "O";
    :param top_number: A value for the top number of de novo sequencing, such as "10";
    :param input_path: A string for the input folder, such as "/data/Training-01-Human-285";
    :return: write all the CSV files to the "/data/Training-01-Human-285/N-GP-CSV";
             write all the denovo msp files to the "/data/Training-01-Human-285/N-GP-DENOVO-MSP".
    """
    path_name = Path(input_path)
    # Check whether path name is a folder
    assert path_name.is_dir(), "Input path is wrong!"

    top_number = int(top_number)
    assert top_number > 0 and top_number < 101, "Wrong Top Number!"

    if input_type == 'N':
        msp_path = path_name / 'N-GP-MSP'
        csv_name = 'N-GP-CSV-TOP-' + str(top_number)
        csv_path = path_name / csv_name
        denovo_name = 'N-GP-DENOVO-MSP-TOP-' + str(top_number)
        denovo_msp_path = path_name / denovo_name
    elif input_type == 'O':
        msp_path = path_name / 'O-GP-MSP'
        csv_name = 'O-GP-CSV-TOP-' + str(top_number)
        csv_path = path_name / csv_name
        denovo_name = 'O-GP-DENOVO-MSP-TOP-' + str(top_number)
        denovo_msp_path = path_name / denovo_name
    else:
        print("Error Input Type!")
        return

    time_start_msp = time.asctime(time.localtime(time.time()))
    start_time = time.time()
    print(f"The start time for parsing msp files: {time_start_msp}")
    print(f"The msp path is: {msp_path.absolute()}")
    print(f"The csv path is: {csv_path.absolute()}")
    print(f"The denovo msp path is: {denovo_msp_path.absolute()}")

    # Make dir for the csv folder
    csv_path.mkdir(parents=True, exist_ok=True)
    denovo_msp_path.mkdir(parents=True, exist_ok=True)

    # Record the total number of files, spectra, spectra for deep learning, and samples.
    total_num_files = 0
    total_num_spectra = 0
    total_num_dl_spectra = 0
    total_num_samples = 0
    # Check the MSP folder
    if not msp_path.exists():
        print("The MSP folder does not exist!")
        return

    # Iterate all the files in the MSP folder
    for msp_file in msp_path.iterdir():
        # Judge whether msp is a folder, only open it as a msp file
        # Example: 202110_Palleon_Cyno_In_vitro-20211027_Palleon_A_HILIC_HCDFT.msp
        if not msp_file.is_file():
            print("There is a folder in the MSP folder!")
            continue
        if msp_file.suffix != '.msp':
           print("There is a file whose type is not msp!")
           continue
        time_read_msp = time.time()
        # Get the filename without the extension.
        msp_stem = msp_file.stem
        # Generate the name for a csv file, and the whole path/file for it.
        csv_name = f'{msp_stem}.csv'
        csv_path_file = csv_path / csv_name
        # Generate the name for a denovo msp file, and the whole path/file for it.
        denovo_name = f'denovo_{msp_stem}.msp'
        denovo_path_file = denovo_msp_path / denovo_name

        num_spectra_file, num_dl_spectra_file, num_samples_file = \
            msp_to_csv_denovo(top_number, msp_file, csv_path_file, denovo_path_file)

        total_num_files += 1
        total_num_spectra += num_spectra_file
        total_num_dl_spectra += num_dl_spectra_file
        total_num_samples += num_samples_file
        time_write_csv = time.time()
        process_time = time_write_csv - time_read_msp
        print(f"The time for reading the msp file, then write the csv file and denovo file: {process_time}")

    time_end_msp =  time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60
    total_time_hours = total_time_minutes / 60
    print(f"The total number of files is: {total_num_files}")
    print(f"The total number of spectra is: {total_num_spectra}")
    print(f"The total number of spectra for the deep learning is: {total_num_dl_spectra}")
    print(f"The total number of samples is: {total_num_samples}")
    print(f"The end time is: {time_end_msp}")
    print(f"The total time is: {total_time_seconds} Seconds, {total_time_minutes} Minutes, {total_time_hours} Hours.")


#  CLI (command line interface) for the input and output
def user_interface(input_type, top_number, input_path):
    """
    Get the folder by the input type and input path for the CSV files
    :param input_type: A string for the folder of input type;
    :param top_number: A value for the top number of de novo sequencing, such as "10";
    :param input_path: A string for the folder of input path.
    :return: The output of CSV files.
    """
    #  python parse_msp_data.py -IT=N -TN=10 -IP=/data/Training-2021-D-Va--Human
    parse_msp_files(input_type, top_number, input_path)


"""
Input parameters for the user interface:
    1   Input Glyco Peptide Type
    2   Input Top Number of De Novo Sequencing
    3   Input File Path Name
"""
parser = argparse.ArgumentParser(description='Input parameters to run the script.')
parser.add_argument('--InputType', '-IT',
                    help='Input Format parameter，required, has default. N-NlinkedGlycoPeptide; O-OlinkedGlycoPeptide.',
                    required=False, default='N')
parser.add_argument('--TopNumber', '-TN',
                    help='Top Number parameter，required, has default. The value of the top number for de novo '
                         'sequencing.', required=False, default='1')
parser.add_argument('--InputPath', '-IP', help='Input Path parameter，required，no default. Such as /data/Training-01-Human-285.',
                    required=False)

args = parser.parse_args()


if __name__ == "__main__":
    # Please run the script with the following input format in Linux/Unix/Mac such as:
    # python parse_msp_to_csv.py -IT=N -TN=10 -IP=/data/Training-01-Human-285
    # Or Windows such as:
    # python parse_msp_data.py -IT=O -TN=10 -IP=D:\\data\\Training-01-Human-285
    # For huge files with several hours, should use "nohup" and "&" to run in the background. For example:
    # nohup python -u parse_msp_to_csv.py -IT=N -TOP=10 -IP=/data/Training-01-Human-285 > Training-01-Human-285_out.out 2>&1 &

    try:
        user_interface(args.InputType, args.TopNumber, args.InputPath)
    except Exception as e:
        print(e)


