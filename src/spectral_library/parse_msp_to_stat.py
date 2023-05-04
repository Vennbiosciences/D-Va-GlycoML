#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read MS/MS spectral library data from an MSP format file

Created on 07 November 2022.

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
    'peptide',
    'peptide_charge',
    'peptide_glycan',
    'peptide_glycan_charge',
]

# Combine 20 amino acids and 5 monosaccharides together to generate 25 codes, the total number is 26 after plus zero.
# An additional amino acid is 'J', which represents the glycosylated cite for 'N'.
# Ignore 'X' which is used to represent both 'I' and 'L'.
# Combine 'I' and 'L' into 'X', save the room for the fifth monosaccharide.
# Using 'I' and 'L' to represent monosaccharides, and 'Z' to represent empty one.
# The original one:  amino_acid_codes = "ACDEFGHIJKLMNPQRSTVWY"
amino_acid_codes = "ACDEFGHJKMNPQRSTVWXY"
zero_code = "Z"
#print(amino_acid_monosaccharide_zero_codes)


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


# read the processed information from a msp file, and write the stat data into a text file.
def msp_to_stat(top_number, msp_file, stat_txt, stat_csv):
    """
    Write the information of input_csv into the files of csv and denovo msp;
    :param top_number: the top number of denovo results;
    :param msp_file: the msp file to read;
    :param stat_txt: the statics txt to write;
    :param stat_csv: the statics csv to write;
    :return: all the statistic information for a msp file.
    """
    num_samples = 0
    num_dl_spectra = 0
    num_spectra = 0

    peptides = set()
    peptide_charges = set()
    peptide_glycans = set()
    peptide_glycan_charges = set()
    glycan_intensity_ratio =[]

    with open(msp_file) as spec_library_file, open(stat_txt, "w", newline='') as stat_output_file, \
            open(stat_csv, "w", newline='') as csv_output_file:
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

                    glycan_intensity = {}
                    total_glcan_intensities = 0.0
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
                            continue
                            # b, or y ions in the peptide backbone
                            # b$2+1, $ means "cross ring fragment of a HexNAc", should be an 83Da
                            # y12-N(1)+2, peptide fragment plus a glycan
                            # For b ions: b3 is 0, b$3 is 1, b3-N(1) is 2,
                            # For y ions: y3 is 3, y$3 is 4, y3-N(1) is 5.
                        elif annotation[0] in glyco_ions:
                            # intact charged peptide without glycan
                            # Y0+2=1113.098
                            if annotation[1] == '0':
                                continue
                            # intact charged peptide with cross ring fragment of a HexNAc
                            # Y$+2=1154.617
                            elif annotation[1] == '$':
                                continue
                            # Y ions: peptide backbone attached a fragmented glycan
                            # Summarize ion intensities for each glycan composition
                            # For GLYCAN(H,N,F,A,G)=5,4,1,0,0 such as
                            # 1018.423	7877.0	Y-H(5)N(4)F(1)+3
                            # annotation[1] == '-':
                            else:
                                glycan_ion = annotation.split("-")[1].split("+")[0]
                                glycan_ion_dic = parse_glycan(glycan_ion)
                                glycan_str = glycan_to_str(glycan_ion_dic)
                                if glycan_str in glycan_intensity:
                                    sum_intensity = round(float(glycan_intensity[glycan_str]) + float(intensity), 1)
                                else:
                                    sum_intensity = round(float(intensity), 1)
                                glycan_intensity[glycan_str] = round(float(sum_intensity), 1)
                                total_glcan_intensities += float(intensity)
                        else:
                            # Ignore the position, use 32 + 18 to replace it
                            # ion_indicator_code = 3
                            # positions.append(50)
                            continue

                    glycan_sequence = glycans.split("=")[1]
                    glycan_dic = parse_glycan(glycan_sequence)
                    glycan_str = glycan_to_str(glycan_dic)
                    # Call de novo sequencing to generate linearized glycan string
                    glycan_deno_lists = denovor.de_novo(glycan_str, glycan_intensity)

                    # Assume the maximum length of glycan is max_glycan_length, pad the left locations with "Z"
                    # The number of spectra for deep learning could have several denovo results for the samples
                    num_dl_spectra += 1
                    # Consider top 10 de novo sequence candidates if the total number greater than 10
                    # top_ten_candidates = min(len(glycan_deno_lists), 10)
                    top_one_candidate = min(len(glycan_deno_lists), top_number)
                    for i in range(top_one_candidate):
                        # glycan_deno_lists[i] is a tuple like ('NNHHNHHANHNHAA', 954586.3)
                        glycan_de_novo_sequence = glycan_deno_lists[i][0]
                        peptide_charge = peptide_name + "_" + charge
                        peptide_glycan = peptide_name + "_" + glycan_de_novo_sequence
                        peptide_glycan_charge = peptide_name + "_" + glycan_de_novo_sequence + "_" + charge
                        glycan_de_novo_intensity = glycan_deno_lists[i][1]
                        peptides.add(peptide_name)
                        peptide_charges.add(peptide_charge)
                        peptide_glycans.add(peptide_glycan)
                        peptide_glycan_charges.add(peptide_glycan_charge)
                        # ignore the zeros
                        if glycan_de_novo_intensity != 0 and total_glcan_intensities != 0:
                            glycan_de_novo_ratio = glycan_de_novo_intensity / total_glcan_intensities
                            glycan_intensity_ratio.append(glycan_de_novo_ratio)
                        else:
                            print("glycan_de_novo_sequence:", glycan_de_novo_sequence)
                            print("total_glcan_intensities:", total_glcan_intensities)
                        num_samples += 1

                    # Should use break here, to successfully end the process of one chunk, and start a new one.
                    break
                else:
                    current_chunk.append(line.strip())
            # If there is no new data, end the reading process.
            if not read:
                break

        glycan_intensity_avg = sum(glycan_intensity_ratio) / len(glycan_intensity_ratio)

        stat_output_file.write(f"The total number of spectra is: {num_spectra} \n")
        stat_output_file.write(f"The total number of spectra for the deep learning is: {num_dl_spectra} \n")
        stat_output_file.write(f"The total number of samples is: {num_samples} \n")
        stat_output_file.write(f"The total number of peptides is: {len(peptides)} \n")
        stat_output_file.write(f"The total number of peptide charges is: {len(peptide_charges)} \n")
        stat_output_file.write(f"The total number of peptide glycans is: {len(peptide_glycans)} \n")
        stat_output_file.write(f"The total number of peptide glycan charges is: {len(peptide_glycan_charges)} \n")
        stat_output_file.write(f"The average intensity ratio of top one denovo to total glycan is: {glycan_intensity_avg} \n")

        # Convert the no duplicated sets to the lists, and align them with equal size,
        peptides_list = list(peptides)
        peptide_charges_list = list(peptide_charges)
        peptide_glycans_list = list(peptide_glycans)
        peptide_glycan_charges_list = list(peptide_glycan_charges)

        # Maximum length: peptide_glycan_charges_list should be the largest one among the four lists
        add_peptides = len(peptide_glycan_charges_list) - len(peptides_list)
        add_peptide_charges = len(peptide_glycan_charges_list) - len(peptide_charges_list)
        add_peptide_glycans = len(peptide_glycan_charges_list) - len(peptide_glycans_list)

        # Add spaces to the end of each list, to make them equal size.
        peptides_list.extend(add_peptides*[''])
        peptide_charges_list.extend(add_peptide_charges*[''])
        peptide_glycans_list.extend(add_peptide_glycans*[''])

        # Iterate four lists, get corresponding values.
        for i in range(len(peptide_glycan_charges_list)):
            row_dictionary = {
                'peptide': peptides_list[i],
                'peptide_charge': peptide_charges_list[i],
                'peptide_glycan': peptide_glycans_list[i],
                'peptide_glycan_charge': peptide_glycan_charges_list[i],
            }
            writer_csv.writerow(row_dictionary)

        #stat_output_file.write(f"The peptides are: \n")
        #stat_output_file.write(str(peptides) + "\n")
        #stat_output_file.write(f"The peptide charges are: \n")
        #stat_output_file.write(str(peptide_charges) + "\n")
        #stat_output_file.write(f"The peptide glycans are: \n")
        #stat_output_file.write(str(peptide_glycans) + "\n")
        #stat_output_file.write(f"The peptide glycan charges are: \n")
        #stat_output_file.write(str(peptide_glycan_charges) + "\n")

        # write the results to the csv file.

    print(f"The total number of spectra is: {num_spectra}")
    print(f"The total number of spectra for the deep learning is: {num_dl_spectra}")
    print(f"The total number of samples is: {num_samples}")
    print(f"Constructed a statistic file: {stat_txt.name}")
    print(f"Constructed a csv file: {stat_csv.name}")
    return num_spectra, num_dl_spectra, num_samples, peptides, peptide_charges, peptide_glycans, \
           peptide_glycan_charges, glycan_intensity_ratio


# Read the path for input folder (MSP) and the top number of de novo sequencing,
# then calculate the total number of peptides, peptide_charges, and peptide_glycan_charges.
def parse_msp_statistic(input_type, top_number, input_path):
    """
    Write the statics file to the folder of "input_path"
    :param input_type: A type for the input folder, such as "N" or "O";
    :param top_number: A value for the top number of de novo sequencing, such as "10";
    :param input_path: A string for the input folder, such as "/data/Training-01-Human-285";
    :return: write a statics file contains the total number of unique peptides and unique peptide-glycansequences, and
    their corresponding sets.
    """
    path_name = Path(input_path)
    # Check whether path name is a folder
    assert path_name.is_dir(), "Input path is wrong!"

    top_number = int(top_number)
    assert top_number > 0 and top_number < 101, "Wrong Top Number!"

    if input_type == 'N':
        msp_path = path_name / 'N-GP-MSP'
        stat_name = f'N-GP-STAT-TOP-{top_number}'
        stat_path = path_name / stat_name
    elif input_type == 'O':
        msp_path = path_name / 'O-GP-MSP'
        stat_name = f'O-GP-STAT-TOP-{top_number}'
        stat_path = path_name / stat_name
    else:
        print("Error Input Type!")
        return

    time_start_msp = time.asctime(time.localtime(time.time()))
    start_time = time.time()
    print(f"The start time for parsing msp files: {time_start_msp}")
    print(f"The msp path is: {msp_path.absolute()}")
    print(f"The stat path is: {stat_path.absolute()}")

    # Make :dir for the statistic folder
    stat_path.mkdir(parents=True, exist_ok=True)

    # Record the total number of files, spectra, spectra for deep learning, and samples.
    total_num_files = 0
    total_num_spectra = 0
    total_num_dl_spectra = 0
    total_num_samples = 0

    # Record the total number of total peptides, total peptide charges, total peptide glycan charges.
    total_peptides = set()
    total_peptide_charges = set()
    total_peptide_glycans = set()
    total_peptide_glycan_charges = set()
    # Record the total intensity ratio of denovo sequence divide glycan.
    total_glycan_intensity_ratio = []

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
        # Generate the name for a stat file, and the whole path/file for it.
        stat_txt = f'{msp_stem}.txt'
        stat_path_txt = stat_path / stat_txt
        stat_csv = f'{msp_stem}.csv'
        stat_path_csv = stat_path / stat_csv

        num_spectra_file, num_dl_spectra_file, num_samples_file, peptides, peptide_charges, peptide_glycans, \
        peptide_glycan_charges, glycan_intensity_ratio = msp_to_stat(top_number, msp_file, stat_path_txt, stat_path_csv)

        total_num_files += 1
        total_num_spectra += num_spectra_file
        total_num_dl_spectra += num_dl_spectra_file
        total_num_samples += num_samples_file

        total_peptides = total_peptides | peptides
        total_peptide_charges = total_peptide_charges | peptide_charges
        total_peptide_glycans = total_peptide_glycans | peptide_glycans
        total_peptide_glycan_charges = total_peptide_glycan_charges | peptide_glycan_charges
        total_glycan_intensity_ratio.extend(glycan_intensity_ratio)

        time_write_stat = time.time()
        process_time = time_write_stat - time_read_msp
        print(f"The time for reading the msp file, then write the statistic file: {process_time}")

    time_end_stat = time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60
    total_time_hours = total_time_minutes / 60
    print(total_glycan_intensity_ratio)
    total_glycan_intensity_avg = sum(total_glycan_intensity_ratio) / len(total_glycan_intensity_ratio)

    print(f"The total number of files is: {total_num_files}")
    print(f"The total number of spectra is: {total_num_spectra}")
    print(f"The total number of spectra for the deep learning is: {total_num_dl_spectra}")
    print(f"The total number of samples is: {total_num_samples}")
    print(f"The total number of peptides is: {len(total_peptides)}")
    print(f"The total number of peptide charges is: {len(total_peptide_charges)}")
    print(f"The total number of peptide glycans is: {len(total_peptide_glycans)}")
    print(f"The total number of peptide glycan charges is: {len(total_peptide_glycan_charges)}")
    print(f"The total average intensity ratio of top one denovo to total glycan is: {total_glycan_intensity_avg} \n")

    #print(f"The peptides are: \n")
    #print(total_peptides)
    #print(f"The peptide charges are: \n")
    #print(total_peptide_charges)
    #print(f"The peptide glycans are: \n")
    #print(total_peptide_glycans)
    #print(f"The peptide glycan charges are: \n")
    #print(total_peptide_g::lycan_charges)

    print(f"The end time is: {time_end_stat}")
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
    parse_msp_statistic(input_type, top_number, input_path)


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


