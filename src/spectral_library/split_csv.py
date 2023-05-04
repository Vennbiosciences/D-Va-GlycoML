#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read a csv file, then split it into multiple csv files with chunk size of 1000000

Created on 24 Aug 2022.
###################################################################################################################
"""
__author__ = 'ZLiang'


import csv
import re
from stack_queue import Stack
from de_novo_sequencing import DeNovoSequencing
import argparse
import pandas as pd
from pathlib import Path
import time


def split_file(input_file, output_path):
    """
    Read a csv file from the input file, then split the csv file into multiple csv files to the output folder,
    :param input_file: A string for the input csv file
    :param output_path: A string for the output csv folder
    :return: write all the processed csv files to output folder.
    """
    time_start_csv = time.asctime(time.localtime(time.time()))
    start_time = time.time()
    print(f"The start time for parsing csv file: {time_start_csv}")

    path_name = Path(output_path)
    # Check whether path name is a folder
    if path_name.is_dir():
        print("Output path is correct!")
    else:
        # Make dir for the output csv folder
        path_name.mkdir(parents=True, exist_ok=True)

    file_num = 0
    for i, chunk in enumerate(pd.read_csv(input_file, chunksize=1000000)):
        output_file = path_name / f"chunk_{i}.csv"
        chunk.to_csv(output_file, index=False)
        file_num += 1
        print(f"Generated a csv file: chunk_{i}.csv")

    time_end_csv =  time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60
    print(f"The total number of files is: {file_num}")
    print(f"The end time is: {time_end_csv}")
    print(f"The total time is: {total_time_seconds} Seconds, {total_time_minutes} Minutes.")


#  CLI (command line interface) for the input and output
def user_interface(input_file, output_path):
    """
    Get the folder by the input type and input path for the CSV files
    :param input_file: A string for the file of input path;
    :param output_path: A string for the folder of input path.
    :return: The output of CSV files.
    """
    #  python parse_msp_data.py -IF=/data/Training-2021-D-Va--Human/N-GP-CSV-TOP-1/GP-ID_HCDFT.csv
    #   -OP=/data/Training-01-Human-285/N-GP-CSV-TOP-1-Split
    split_file(input_file, output_path)


"""
Input parameters for the user interface:
    1   Input File Path Name
    2   Output Path Name
"""
parser = argparse.ArgumentParser(description='Input parameters to run the script.')
parser.add_argument('--InputFile', '-IF', help='Input File parameter，required，no default. Such as '
                                    '/data/Training-01-Human-285/N-GP-CSV-TOP-1/GP-ID_HCDFT.csv.', required=False)
parser.add_argument('--OutputPath', '-OP', help='Output Path parameter，required，no default. Such as '
                                    '/data/Training-01-Human-285/N-GP-CSV-TOP-1-Split.', required=False)
args = parser.parse_args()


if __name__ == "__main__":
    # Please run the script with the following input format in Linux/Unix/Mac such as:
    # python split_csv.py -IF=/data/Training-01-Human-285/N-GP-CSV-TOP-1/GP-ID_HCDFT.csv
    # -OP=/data/Training-01-Human-285/N-GP-CSV-TOP-1-Split
    # Or Windows such as:
    # python split_csv.py -IP=D:\\data\\Training-01-Human-285\\N-GP-CSV-TOP-1\\GP-ID_HCDFT.csv
    # -OP=D:\\data\\Training-01-Human-285\\N-GP-CSV-TOP-1-Split

    try:
        user_interface(args.InputFile, args.OutputPath)
    except Exception as e:
        print(e)

