#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
#######################################################################################################################
This script read all the original csv files in a folder, whose path is from user's input.
Then converts each row into a data format for X, X_meta, and y for the corresponding deep learning model training.
These data sets could be stored in two kinds of file formats: csv and pkl.
For a deep learning csv file, which is readable for the human, the size is triple than that of an original csv file,
and it will spend one to two minutes to finish writing.
The reason is for the empty spots of b, y, and Y ions, still need pad 0s for their intensities.
For the pkl, which converts the information into a byte stream, and writes it into a pickle file, that implements
binary protocols for serializing.
For the next step of training and testing, unpickling is the inverse operation, whereby a binary file is converted
back into a format for X, X_meta, and y.
The size of a pkl file is about triple than that of a csv file, the speed for writing is only one second.
But the speed of reading pkl file is faster than that of csv file, roughly double.

Created on 8 March 2022.
Modified on 11 April 2022, to handle both N and O linked glycopeptides.
#######################################################################################################################
"""
__author__ = 'ZLiang'

import csv
import numpy as np
import pandas as pd
import pickle
from pathlib import Path
import argparse
import time

# Control whether to truncate or not. Should set as infinity "inf".
# Otherwise, it outputs ellipsis "..." to represent the truncated lists, which affects csv files, but not for pkl files.
np.set_printoptions(threshold=np.inf)

# Assume the maximum length of peptide is 32, and maximum length of glycan is 18.
# The maximum charge of glycopeptide is 6.
# The maximum number of ions is 3 (b, y, Y).
max_peptide_length = 32
max_glycan_length = 18
max_num_charges = 6
max_num_ions = 3

# sequence length
seq_len = max_peptide_length + max_glycan_length

# Header for the deep learning input output csv
dl_csv_header = [
    'X, X_meta, y',
]

# write the processed information into a pickle file, should use binary method: 'b'.
def data_to_pkl(input_pkl, output_pkl_file):
    """
    Write the information of input_pkl into the file of out_pkl_file
    :param input_pkl: the data of X, X_meta, y.
    :param output_pkl_file: the path and name of the pickle file.
    :return:
    """

    with open(output_pkl_file, 'wb') as pkl_file:
        pickle.dump(input_pkl, pkl_file)
    pkl_name = output_pkl_file.name
    print(f"Constructed a pkl file:{pkl_name}")


# write the processed information into a deep learning input and output csv file.
def data_to_csv(input_csv, output_dl_csv_file):
    """
    Write the information of input_csv into the file of output_dl_csv_file
    :param input_csv: the data of X, X_meta, y;
    :param output_dl_csv_file: the path and name of the deep learning csv file.
    :return:
    """

    with open(output_dl_csv_file, "w", newline='') as dl_csv_file:
        writer = csv.DictWriter(dl_csv_file, fieldnames=dl_csv_header)
        writer.writeheader()
        for line in input_csv:
            row_dictionary = {
                 'X, X_meta, y': line,
            }
            writer.writerow(row_dictionary)

    dl_csv_name = output_dl_csv_file.name
    print(f"Constructed a dl csv file:{dl_csv_name}")


# read the processed information from a csv file, and write the data into a pkl file and a deep learning csv file.
def csv_to_pkl_dl(csv_file, pkl_file, dl_csv_file):
    """
    Write the information of input_csv into the files of csv and denovo msp;
    :param csv_file: the csv file to read;
    :param pkl_file: the pkl file to write;
    :param dl_csv_file: the deep leaning csv file to write.
    :return: the total number of samples for a pkl file.
    """
    num_samples = 0
    df = pd.read_csv(csv_file, encoding='utf-8')
    # Initialize X, y, X_meta for the deep learning model.
    Xs = []
    ys = []
    X_metas = []
    for line in np.array(df):
        # Data for X: glycopeptide one-hot encoding (50, 26), use "eval" to keep it as lists.
        Xs.append(np.array(eval(line[0])))
        # Data for X_meta: charge (1), should use np.array to generate an array, add [] for metas.
        X_metas.append(np.array([line[1]]))
        # Data for y: (50, 18)
        # intensities
        intensities = np.array(eval(line[5])).astype(float)
        # ions
        ions = np.array(eval(line[6])).astype(int)
        # positions
        positions = np.array(eval(line[7])).astype(int)
        # ion charges
        ion_charges = np.array(eval(line[8])).astype(int)

        assert np.min(intensities) > 0
        assert np.min(ions) >= 0
        assert np.min(positions) >= 1
        assert np.min(ion_charges) >= 0

        # Consider four different types of ions: b ion, y ion, Y ion.
        # The max position of the ion is limited to "seq_len - 1".
        # For the length glycopeptide of 50, we have 49 for y.
        # Use "b, y, Y" to represent ions, and times the maximum charges, then the space is 3*6 = 18
        y = np.zeros((seq_len - 1, max_num_ions * max_num_charges))
        # Arrange ion intensities to the corresponding places.
        for i in range(len(intensities)):
            # b ion, keep the positions[i], and use [1, 0, 0] for the b ion.
            if ions[i] == '0':
                position = positions[i] - 1
                ion_type = (ion_charges[i] - 1) * max_num_ions + 0
            # y ion, use max_peptide_length - positions[i], and use [0, 1, 0] for the y ion.
            elif ions[i] == '1':
                position = max_peptide_length - positions[i] - 1
                ion_type = (ion_charges[i] - 1) * max_num_ions + 1
            # Y ion, ions[i] == '2', keep the positions[i], and use [0, 0, 1] for the Y ion.
            else:
                position = positions[i] - 1
                ion_type = (ion_charges[i] - 1) * max_num_ions + 2
            # Add the intensity to the corresponding place.
            y[position, ion_type] += intensities[i]
        # Normalization by the summary of the total intensities.
        y = y / np.sum(y)
        # Concatenate into y(50, 18)
        ys.append(y)
        num_samples += 1

    # Construct the final data for the deep learning model.
    dl_data = list(zip(Xs, X_metas, ys))

    time_write_pkl = time.asctime(time.localtime(time.time()))
    print(f"The time for writing the pkl file:{time_write_pkl}")
    # write the data into a pickle file.
    data_to_pkl(dl_data, pkl_file)

    time_write_dl_csv = time.asctime(time.localtime(time.time()))
    print(f"The time for writing the dl csv file:{time_write_dl_csv}")
    # write the data into a deep learning input and output csv file.
    data_to_csv(dl_data, dl_csv_file)

    return num_samples


# Read a path for the input folder (CSV), then process the csv files into X, X_meta, and y,
# which are used in the deep learning model.
def parse_csv_files(input_type, input_path):
    """
    Read the CSV files from the folder of "input_path", then converts them into the formats for the deep learning model
    :param input_type: A string for the folder of input type, such as "N" or "O"
    :param input_path: A string for the input folder, such as "/data/Training-01-Human-285".
    :return: write all the PKL files to the "/data/Training-01-Human-285/N-GP-PKL";
             write all the denovo msp files to the "/data/Training-01-Human-285/N-GP-DENOVO-MSP".
    """
    path_name = Path(input_path)

    # Check whether path name is a folder
    assert path_name.is_dir(), "Input path is wrong!"

    if input_type == 'N':
        csv_path = path_name / 'N-GP-CSV'
        pkl_path = path_name / 'N-GP-PKL'
        dl_csv_path = path_name / 'N-GP-DL-CSV'
    elif input_type == 'O':
        csv_path = path_name / 'O-GP-CSV'
        pkl_path = path_name / 'O-GP-PKL'
        dl_csv_path = path_name / 'O-GP-DL-CSV'
    else:
        print("Error Input Type!")
        return

    time_start_csv = time.asctime(time.localtime(time.time()))
    start_time = time.time()
    print(f"The start time for parsing csv files: {time_start_csv}")
    print(f"The csv path is: {csv_path.absolute()}")
    print(f"The pkl path is: {pkl_path.absolute()}")
    print(f"The deep learning csv path is: {dl_csv_path.absolute()}")

    # Make dir for the pickle folder
    pkl_path.mkdir(parents=True, exist_ok=True)
    # Make dir for the deep learning csv folder
    dl_csv_path.mkdir(parents=True, exist_ok=True)

    # Record the total number of files and rows.
    total_num_files = 0
    total_num_samples = 0

    # Explore all the files in the CSV folder
    if not csv_path.exists():
        print("The CSV folder does not exist!")
        return

    # Iterate all the files in the CSV folder
    for csv_file in csv_path.iterdir():
        if not csv_file.is_file():
            print("There is a folder in the CSV folder!")
            continue
        if csv_file.suffix != '.csv':
           print("There is a file whose type is not csv!")
           continue
        time_read_csv = time.asctime(time.localtime(time.time()))
        print(f"The time for reading the csv file: {time_read_csv}")

        # Get the filename without the extension, csv_file is the full "path/name" for the csv file.
        csv_stem = csv_file.stem
        # Generate the name for a pkl file, and the whole path/file for it.
        pkl_name = f'{csv_stem}.pkl'
        pkl_path_file = pkl_path / pkl_name
        # Generate the name for deep learning input and output csv file, and the whole path/file for it.
        dl_csv_name = f'dl_{csv_stem}.csv'
        dl_csv_path_file = dl_csv_path / dl_csv_name

        num_samples_one_file = csv_to_pkl_dl(csv_file, pkl_path_file, dl_csv_path_file)

        total_num_files += 1
        total_num_samples += num_samples_one_file

    time_end_csv = time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60
    total_time_hours = total_time_minutes / 60
    print(f"The total number of files is: {total_num_files}")
    print(f"The total number of samples is: {total_num_samples}")
    print(f"The end time is: {time_end_csv}")
    print(f"The total time is: {total_time_seconds} Seconds, {total_time_minutes} Minutes, {total_time_hours} Hours.")

#  CLI (command line interface) for the input and output
def user_interface(input_type, input_path):
    """
    Get the folder of input path for the CSV and PKL files
    :param input_type: A string for the folder of input type;
    :param input_path: A string for the folder of input path.
    :return: Call read_csv for the input of CSV files, and convert them into PKL and deep learning CSV files.
    """
    #  python parse_csv_data.py -IP=/data/Training-2021-D-Va--Human
    parse_csv_files(input_type, input_path)

"""
Input parameters for the user interface:
    1   Input Glyco Peptide Type
    2   Input File Path Name
"""
parser = argparse.ArgumentParser(description='Input parameters to run the script.')
parser.add_argument('--InputType', '-IT',
                    help='Input Format parameter，required, has default. N-NlinkedGlycoPeptide; O-OlinkedGlycoPeptide.',
                    required=False, default='N')
parser.add_argument('--InputPath', '-IP',
                    help='Input Path parameter，required，no default. Such as /data/Training-01-Human-285.',
                    required=False)

args = parser.parse_args()

if __name__ == "__main__":
    # Please run the script with the following input format in Linux/Unix/Mac such as:
    # python parse_csv_data.py -IT=N -IP=/data/Training-01-Human-285
    # Or Windows such as:
    # python parse_csv_data.py -IT=O -IP=D:\data\Training-01-Human-285
    # For huge files processed by several hours, should use "nohup" and "&" to run in the background. For example:
    # nohup python -u parse_csv_data.py -IP=/data/Training-2021-D-Va--Human > D-Va–Human_out_pkl.out 2>&1 &

    try:
        user_interface(args.InputType, args.InputPath)
    except Exception as e:
        print(e)




