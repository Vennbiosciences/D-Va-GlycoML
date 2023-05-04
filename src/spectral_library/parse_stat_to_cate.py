#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read a statistic csv file from, then divide them into three N-glycoform categorizations.

Info from Matthew Campbell:
If the number of HexNAc residues equals 2 and the number of hexose residues is greater than or equal to 5,
    then the N-linked glycan is of the type “high mannose”;
If the number of HexNAc residues is greater than or equal to 3 and the number of hexose residues is also greater than
    or equal to 3, then the N-linked glycan is of the type “hybrid/complex”;
If the number of HexNAc residues is equal to 2 and the number of hexose residue residues is less than or equal to 3,
    then the N-linked glycan is of the type “paucimannose”.

Created on 14 November 2022.

###################################################################################################################
"""
__author__ = 'ZLiang'

import csv
import argparse
from pathlib import Path
import time
import pandas as pd

monosaccharide_list = "HNFAG"

# Assume the maximum length of peptide is 32, and maximum length of glycan is 18.
MAX_PEPTIDE_LENGTH = 32
MAX_GLYCAN_LENGTH = 18

# Using characters to represent five monosaccharides
character_monosaccharide_replacements = {
    "H": "Hex",
    "N": "HexNAc",
    "F": "Fuc",
    "A": "NeuAc",
    "G": "NeuGc"
}

# Headers for the csv output
csv_output_rows = [
    'peptide_glycan',
    'peptide_glycan_categorization',
    'peptide_glycan_charge',
    'peptide_glycan_charge_categorization',
]


def sequence_to_list(glycan_seq):
    """
    Consider double-digit monosaccharide, use two digits to represent a monosaccharide.
    Parse glycan sequence such as "HHNNF" into a list [2, 2, 1, 0, 0], which represents the number for "HNFAG"
    :param glycan_dic: A sting like {character+}
    :return: A number list like [2, 2, 1, 0, 0]
    """
    glycan_mono = "HNFAG"
    glycan_number = [0, 0, 0, 0, 0]
    for char in glycan_seq:
        glycan_number[glycan_mono.find(char)] += 1

    return glycan_number

def glycan_to_categorization(glycan_number):
    """
    Consider double-digit monosaccharide, use two digits to represent a monosaccharide.
    Parse glycan number list such as [2, 2, 1, 0, 0], into three N-linked glycan categorizations
    :param glycan_seq: A number list like [2, 2, 1, 0, 0]
    :return: A categorization from {'', 'High Mannose', 'Hybrid/Complex', 'Paucimannose'}
    """
    # glycan_number is a number list for 'HNFAG'
    # if 'N' == 2 and 'H' >= 5, then it is “high mannose” N-linked glycan
    if glycan_number[1] == 2 and glycan_number[0] >= 5:
        return 'High Mannose'
    # if 'N' >= 3 and 'H' >= 3, then it is “hybrid/complex” N-linked glycan
    elif glycan_number[1] >= 3 and glycan_number[0] >= 3:
        return 'Hybrid/Complex'
    # if 'N' ==  and 'H' <= 3, then it is “hybrid/complex” N-linked glycan
    elif glycan_number[1] == 2 and glycan_number[0] <= 3:
        return 'Paucimannose'
    else:
        return ''

# read the processed information from a statistic csv file, and write the categorization data into a text file.
def stat_to_cate(stat_csv, cate_txt, cate_csv):
    """
    Write the information of input_csv into the files of csv and denovo msp;
    :param top_number: the top number of denovo results;
    :param stat_csv: the csv file to read;
    :param cate_txt: the categorization txt to write;
    :param cate_csv: the categorization csv to write;
    :return: all the categorization information for a statistic csv file.
    """
    num_glycopeptides = 0
    num_glycopeptides_high_mannose = 0
    num_glycopeptides_hybrid_complex = 0
    num_glycopeptides_paucimannose = 0
    num_glycopeptide_charges = 0
    num_glycopeptide_charges_high_mannose = 0
    num_glycopeptide_charges_hybrid_complex = 0
    num_glycopeptide_charges_paucimannose = 0

    with open(stat_csv) as stat_input_file, open(cate_txt, "w", newline='') as cate_output_file, open(cate_csv, "w", newline='') as csv_output_file:
        df_stat = pd.read_csv(stat_input_file, low_memory=False, sep=',', header=0, encoding='utf8')
        writer_csv = csv.DictWriter(csv_output_file, fieldnames=csv_output_rows)
        writer_csv.writeheader()

        for index, row in df_stat.iterrows():
            # GGSSGWSGGLAQJR_NNHHHNHNHAA
            peptide_glycan = str(df_stat.iloc[index]['peptide_glycan'])
            # After read as dataframe, '' in csv will be converted as 'nan'
            if peptide_glycan != 'nan':
                glycan_sequence = peptide_glycan.split('_')[1]
                glycan_number = sequence_to_list(glycan_sequence)
                glycan_categorization = glycan_to_categorization(glycan_number)
                num_glycopeptides += 1
                if glycan_categorization == 'High Mannose':
                    num_glycopeptides_high_mannose += 1
                elif glycan_categorization == 'Hybrid/Complex':
                    num_glycopeptides_hybrid_complex += 1
                elif glycan_categorization == 'Paucimannose':
                    num_glycopeptides_paucimannose += 1
            else:
                peptide_glycan = ''
                glycan_categorization = ''
            # ADGTVNQIEGEATPVJLTEPAK_NNHHHFNAHHN_4
            peptide_glycan_charge = str(df_stat.iloc[index]['peptide_glycan_charge'])
            glycan_charge_sequence = peptide_glycan_charge.split('_')[1]
            glycan_charge_number = sequence_to_list(glycan_charge_sequence)
            glycan_charge_categorization = glycan_to_categorization(glycan_charge_number)
            num_glycopeptide_charges += 1
            if glycan_charge_categorization == 'High Mannose':
                num_glycopeptide_charges_high_mannose += 1
            elif glycan_charge_categorization == 'Hybrid/Complex':
                num_glycopeptide_charges_hybrid_complex += 1
            elif glycan_charge_categorization == 'Paucimannose':
                num_glycopeptide_charges_paucimannose += 1

            row_dictionary = {
                'peptide_glycan': peptide_glycan,
                'peptide_glycan_categorization': glycan_categorization,
                'peptide_glycan_charge': peptide_glycan_charge,
                'peptide_glycan_charge_categorization': glycan_charge_categorization,
            }
            writer_csv.writerow(row_dictionary)

        cate_output_file.write(f"The number of glycopeptides is: {num_glycopeptides} \n")
        cate_output_file.write(f"The number of glycopeptides high mannose is: {num_glycopeptides_high_mannose} \n")
        cate_output_file.write(f"The number of glycopeptides hybrid complex is: {num_glycopeptides_hybrid_complex} \n")
        cate_output_file.write(f"The number of glycopeptides paucimannose is: {num_glycopeptides_paucimannose} \n")
        cate_output_file.write(f"The number of glycopeptide charges is: {num_glycopeptide_charges} \n")
        cate_output_file.write(f"The number of glycopeptide charges high mannose is: {num_glycopeptide_charges_high_mannose} \n")
        cate_output_file.write(f"The number of glycopeptide charges hybrid complex is: {num_glycopeptide_charges_hybrid_complex} \n")
        cate_output_file.write(f"The number of glycopeptide charges paucimannose is: {num_glycopeptide_charges_paucimannose} \n")

    print(f"The number of glycopeptides is: {num_glycopeptides} \n")
    print(f"The number of glycopeptides high mannose is: {num_glycopeptides_high_mannose} \n")
    print(f"The number of glycopeptides hybrid complex is: {num_glycopeptides_hybrid_complex} \n")
    print(f"The number of glycopeptides paucimannose is: {num_glycopeptides_paucimannose} \n")
    print(f"The number of glycopeptide charges is: {num_glycopeptide_charges} \n")
    print(f"The number of glycopeptide charges high mannose is: {num_glycopeptide_charges_high_mannose} \n")
    print(f"The number of glycopeptide charges hybrid complex is: {num_glycopeptide_charges_hybrid_complex} \n")
    print(f"The number of glycopeptide charges paucimannose is: {num_glycopeptide_charges_paucimannose} \n")

    return num_glycopeptides, num_glycopeptide_charges


# Read the path for input folder (CSV) and the top number of de novo sequencing,
# then calculate three categorizations of glycopeptides.
def parse_statistic_categorization(input_type, top_number, input_path):
    """
    Write the statics file to the folder of "input_path"
    :param input_type: A type for the input folder, such as "N" or "O";
    :param top_number: A value for the top number of de novo sequencing, such as "10";
    :param input_path: A string for the input folder, such as "/data/Training-01-Human-285";
    :return: write a categorization file for each of the glycopeptide.
    """
    path_name = Path(input_path)
    # Check whether path name is a folder
    assert path_name.is_dir(), "Input path is wrong!"

    top_number = int(top_number)
    assert top_number > 0 and top_number < 101, "Wrong Top Number!"

    if input_type == 'N':
        stat_name = f'N-GP-STAT-TOP-{top_number}'
        stat_path = path_name / stat_name
        cate_name = f'N-GP-CATE-TOP-{top_number}'
        cate_path = path_name / cate_name
    elif input_type == 'O':
        stat_name = f'O-GP-STAT-TOP-{top_number}'
        stat_path = path_name / stat_name
        cate_name = f'O-GP-CATE-TOP-{top_number}'
        cate_path = path_name / cate_name
    else:
        print("Error Input Type!")
        return

    time_start_stat = time.asctime(time.localtime(time.time()))
    start_time = time.time()
    print(f"The start time for parsing statistic files: {time_start_stat}")
    print(f"The statistic path is: {stat_path.absolute()}")
    print(f"The categorization path is: {cate_path.absolute()}")

    # Make :dir for the statistic folder
    cate_path.mkdir(parents=True, exist_ok=True)

    # Record the total number of files, glyco peptide ID.
    total_num_files = 0
    total_num_glycopeptide_ids = 0
    total_num_glycopeptide_charge_ids = 0

    # Check the MSP folder
    if not cate_path.exists():
        print("The statistic folder does not exist!")
        return

    # Iterate all the files in the MSP folder
    for stat_file in stat_path.iterdir():
        # Judge whether msp is a folder, only open it as a csv file
        # Example: 202110_Palleon_Cyno_In_vitro-20211027_Palleon_A_HILIC_HCDFT.csv
        if not stat_file.is_file():
            print("There is a folder in the statistic folder!")
            continue
        if stat_file.suffix != '.csv':
            print("There is a statistic text file whose!")
            continue
        time_read_msp = time.time()
        # Get the filename without the extension.
        stat_stem = stat_file.stem
        # Generate the name for a cate file, and the whole path/file for it.
        cate_txt = f'{stat_stem}.txt'
        cate_csv = f'{stat_stem}.csv'
        cate_path_txt = cate_path / cate_txt
        cate_path_csv = cate_path / cate_csv

        num_glycopeptide_ids, num_glycopeptide_charge_ids = stat_to_cate(stat_file, cate_path_txt, cate_path_csv)

        total_num_files += 1
        total_num_glycopeptide_ids += num_glycopeptide_ids
        total_num_glycopeptide_charge_ids += num_glycopeptide_charge_ids

        time_write_stat = time.time()
        process_time = time_write_stat - time_read_msp
        print(f"The time for reading the msp file, then write the statistic file: {process_time}")

    time_end_cate = time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60
    total_time_hours = total_time_minutes / 60
    print(f"The total number of files is: {total_num_files}")
    print(f"The total number of glycopeptide ids is: {total_num_glycopeptide_ids}")
    print(f"The total number of glycopeptide charge ids is: {total_num_glycopeptide_charge_ids}")

    print(f"The end time is: {time_end_cate}")
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
    parse_statistic_categorization(input_type, top_number, input_path)


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
parser.add_argument('--InputPath', '-IP',
                    help='Input Path parameter，required，no default. Such as /data/Training-01-Human-285.',
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
