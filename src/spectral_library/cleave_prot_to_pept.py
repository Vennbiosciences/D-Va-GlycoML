#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read protein sequences from a protein fasta database file,
then cleaves each protein with trypsin digestion, and output all the corresponding peptides.

Created on 24 January 2023.

###################################################################################################################
"""
__author__ = 'ZLiang'


import csv
import argparse
from pathlib import Path
from pyteomics.parser import xcleave
from pyteomics import fasta

# Headers for the csv output
csv_output_rows = [
    'peptide',
    'protein_id',
    'protein_name',
    'index',
]

# Read the fasta file as the input, then cleave each protein with trypsin digestion,
# and write all the digested peptides into a csv file.
def cleave_protein_peptides(input_database_file, output_peptide_path):
    """
    Write a CSV file with all digested peptides from the dababase.
    :param input_database_file: A string for the folder of input fasta file.
    :param output_peptide_path: A string for the output csv folder
    :return: write all the digested peptides into a specified csv file.
    """
    database_name = Path(input_database_file)
    # Check whether database name is a file
    assert database_name.is_file(), "Input database file is wrong!"

    database_stem = database_name.stem
    peptide_csv = f'{database_stem}.csv'

    peptide_path = Path(output_peptide_path)
    # Make dir for the output folder
    peptide_path.mkdir(parents=True, exist_ok=True)

    # Get the path for the output csv file
    peptide_path_csv = peptide_path / peptide_csv

    # lists for the fasta information
    fasta_id_list = []
    fasta_prot_list = []
    fasta_seq_list = []

    # total number of peptides
    num_peptides = 0

    with open(database_name) as database_file, open(peptide_path_csv, "w", newline='') as csv_output_file:
        writer_csv = csv.DictWriter(csv_output_file, fieldnames=csv_output_rows)
        writer_csv.writeheader()

        for description, sequence in fasta.FASTA(database_file):
            prot_id = description.split("|")[1]
            prot_name = description.split("|")[2].split(" ")[0]
            fasta_id_list.append(prot_id)
            fasta_prot_list.append(prot_name)
            fasta_seq_list.append(sequence)

            # parameters for cleave:
            # missed_cleavages: 1
            # min_length: Minimum peptide length
            # max_length: Maximum peptide length
            peptide_list = xcleave(sequence, 'trypsin', 1, 6, 32)

            for peptide in peptide_list:
                row_dictionary = {
                    'peptide': peptide[1],
                    'protein_id': prot_id,
                    'protein_name': prot_name,
                    'index': peptide[0],
                }
                writer_csv.writerow(row_dictionary)
                num_peptides += 1

    print(f"The total number of peptides is: {num_peptides}")
    return num_peptides


#  CLI (command line interface) for the input and output
def user_interface(input_database_file, output_peptide_path):
    """
    Get the fasta file as the input, and output the cleavages into the CSV files
    :param input_database_file: A string for the path of input fasta file.
    :param output_peptide_path: A string for the folder of output csv file.
    :return: The output of CSV files.
    """
    #  python cleave_prot_to_pept.py -IDF=/data/MS-DL-04-D-Va-modeling-22-07-v1.1/database/Human-H.sapiens-SP-1808.fasta
    #      -OPF=/data/MS-DL-04-D-Va-modeling-22-07-v1.1/DigestedPeptides/Human-H.sapiens-SP-1808.csv
    cleave_protein_peptides(input_database_file, output_peptide_path)


"""
Input parameters for the user interface:
    1   Input Database File
    2   Output Peptide Path
"""
parser = argparse.ArgumentParser(description='Input parameters to run the script.')
parser.add_argument('--InputDatabaseFile', '-IDF', help='Input database file parameter，required，no default. '
                        'Such as /data/MS-DL-04-D-Va-modeling-22-07-v1.1/Database/Human-H.sapiens-SP-1808.fasta.',
                    required=False)
parser.add_argument('--OutputPeptidePath', '-OPP', help='Output digested peptide path parameter，required，no default. '
                        'Such as /data/MS-DL-04-D-Va-modeling-22-07-v1.1/DigestedPeptides/.',
                    required=False)
args = parser.parse_args()


if __name__ == "__main__":
    # Please run the script with the following input format in Linux/Unix/Mac such as:
    # python cleave_prot_to_pept.py -IDF=/data/MS-DL-04-D-Va-modeling-22-07-v1.1/database/Human-H.sapiens-SP-1808.fasta
    #   -OPP=/data/MS-DL-04-D-Va-modeling-22-07-v1.1/DigestedPeptides/
    # Or Windows such as:
    # python cleave_prot_to_pept.py -IDF=D:\\data\\database\\Human-H.sapiens-SP-1808.fasta\\
    #   -OPP=D:\\data\\DigestedPeptides\\

    try:
        user_interface(args.InputDatabaseFile, args.OutputPeptidePath)
    except Exception as e:
        print(e)