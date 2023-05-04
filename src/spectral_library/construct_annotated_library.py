#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
################################################################################
This script combines mgf, pGlyco3, and gLabel to generate annotated glycopeptide spectral libraries.
The glycopeptide spectral library could be used for machine learning and deep learning.

For peptide backbone, keep all the b-ions and y-ions.
For intact glycopeptide, only keep Y-ions.

Created on 10 August 2021
Modified On 15 Feb 2022 for five monosaccharides: glycan(H,N,F,A,G), add G here.
Add Command Line Interface (CLI) On 16 Feb 2022 for the user input.
Modified On 27 Mar 2022: add Glycopeptide for "N" or "O" in the CLI.
Modified on 13 April 2022, for the handling of N and O linked glycopeptides.
Modified on 22 Aug 2022, for the handling of duplicated IDs from gLabel that matched a spectrum from mgf file.
################################################################################
"""
__author__ = 'Zhewei Liang: zliang@venn.bio'


from pyteomics import mgf
import pandas as pd
import argparse
from pathlib import Path
import time


# Binary search for the glabel['spec']
def binary_search(spectra, target):
    """
    using binary search to find the target from the spectra
    :param spectra:
    :param target:
    :return:
    """
    if len(spectra) == 0:
        return -1
    left = mid = 0
    right = len(spectra) - 1
    while left <= right:
        # should use int here
        mid = int(left + (right - left) / 2)
        if spectra.iloc[mid] == target:
            return mid
        elif spectra.iloc[mid] < target:
            left = mid + 1
        else:
            right = mid - 1
    # End condition: left > right
    return -1


# Read the information from a mgf and gLabel file, then write the spectral library into a msp file.
def mgf_glabel_to_msp(mgf_file, glabel_file, msp_file):
    """
    Write the spectral library data from mgf and glabel, into a msp file;
    :param mgf_file: the mgf file to read;
    :param glabel_file: the glabel file to read;
    :param msp_file: the msp file to write.
    :return: the total number of samples for a mgf file.
    """
    num_samples = 0
    with open(msp_file, 'w') as msp_writer:
        # Load gLabel file
        glabel = pd.read_csv(glabel_file, sep='\t', header=0, encoding='utf8', engine='python', on_bad_lines='skip')
        glabel_sort = glabel.sort_values(by=["spec"])
        # Need to reset index after sorting
        glabel_sort = glabel_sort.reset_index(drop=True)
        # mgf is a Path Object, should use "str" to convert it for the pyteomics.
        with mgf.read(str(mgf_file.absolute())) as mgf_reader:
            for spectrum in mgf_reader:
                index = binary_search(glabel_sort['spec'], spectrum['params']['title'])
                # Could not find it
                if index == -1:
                    continue
                # One spectrum from MGF might have several identifications in gLabel.
                # After find the first matching, we could continue to search the others in the following orders.
                first = True
                while spectrum['params']['title'] == glabel_sort['spec'][index]:
                    # Generate an entry for the spectral library
                    anno_array = ['' for _ in range(len(spectrum['m/z array']))]
                    matched_ions = glabel_sort.iloc[index]['matched_ion'].split(';')
                    ion_dic = {}
                    for ion in matched_ions:
                        # ion format is 'Y-H(3)N(3)+1=1996.916,4248.5', ion_anno[0] is Y-H(3)N(3)+1, ion_mass[0] is 1996.916
                        ion_anno = ion.split('=')
                        ion_mass = ion_anno[1].split(',')
                        ion_dic[ion_mass[0]] = ion_anno[0]
                    # for each MS2, find the annotation information，
                    for i in range(len(spectrum['m/z array'])):
                        mz = spectrum['m/z array'][i]
                        mz_str = str(mz)
                        if mz_str in ion_dic:
                            anno_array[i] = ion_dic.get(mz_str)
                    # For the MS2 spectra,
                    # if annotated, keep it for the new MS2 information, otherwise remove it now
                    spectrum_mz = []
                    spectrum_intensity = []
                    spectrum_annotation = []
                    for i in range(len(spectrum['m/z array'])):
                        if anno_array[i] != '':
                            spectrum_mz.append(spectrum['m/z array'][i])
                            spectrum_intensity.append(spectrum['intensity array'][i])
                            spectrum_annotation.append(anno_array[i])
                    msp_writer.write('Name: ' + spectrum['params']['title'] + '\n')
                    msp_writer.write('MW: ' + str(spectrum['params']['pepmass'][0]) + '\n')
                    msp_writer.write('Comment: ' + 'CHARGE=' + str(spectrum['params']['charge']) + '\t' +
                                     'RTINSECONDS=' + str(spectrum['params']['rtinseconds']) + '\t' +
                                     'PEPTIDE=' + glabel_sort.iloc[index]['peptide'] + '\t' +
                                     'MODIFICATIONS=' + str(glabel_sort.iloc[index]['modinfo']) + '\t' +
                                     'GLYCAN(H,N,F,A,G)=' + str(glabel_sort.iloc[index]['glycan(H,N,F,A,G)']) + '\t' +
                                     'GLYCANS=' + str(glabel_sort.iloc[index]['formula']) + '\n')
                    msp_writer.write(f'Num Peaks: {len(spectrum_mz)}' + '\n')
                    for i in range(len(spectrum_mz)):
                        msp_writer.write(str(spectrum_mz[i]) + '\t' + str(spectrum_intensity[i]) + '\t'
                                         + str(spectrum_annotation[i]) + '\n')
                    msp_writer.write('\n')
                    num_samples += 1
                    if glabel_sort['spec'][index] == glabel_sort['spec'][index + 1]:
                        print("Found duplicates:", glabel_sort['spec'][index])
                    index += 1
        print(f"Constructed a msp file: {msp_file}")
    return num_samples


# Read the path for input folders (MGF, and gLabel), then write the msp files into MSP folder.
def parse_mgf_files(input_type, input_path):
    """
    Write the MSP files to the folder of "input_path"
    :param input_type: A type for the input folder, such as "N" or "O"
    :param input_path: A string for the input folder, such as "/data/Training-01-Human-285".
    :return: write all the MSP files to the "/data/Training-01-Human-285/N-GP-MSP".
    """
    path_name = Path(input_path)
    # Check whether path name is a folder
    assert path_name.is_dir(), "Input path is wrong!"

    mgf_path = path_name / 'MGF'
    if input_type == 'N':
        glabel_path = path_name / 'N-GP-gLabel'
        msp_path = path_name / 'N-GP-MSP'
    elif input_type == 'O':
        glabel_path = path_name / 'O-GP-gLabel'
        msp_path = path_name / 'O-GP-MSP'
    else:
        print("Error Input Type!")
        return

    time_start_mgf = time.asctime(time.localtime(time.time()))
    start_time = time.time()
    print(f"The start time for parsing mgf files: {time_start_mgf}")
    print(f"The mgf path is: {mgf_path.absolute()}")
    print(f"The gLabel path is: {glabel_path.absolute()}")
    print(f"The msp path is: {msp_path.absolute()}")

    # Make dir for the msp folder
    msp_path.mkdir(parents=True, exist_ok=True)

    # Record the total number of files and samples.
    total_num_files = 0
    total_num_samples = 0
    # Check the MGF and gLabel folder
    if not mgf_path.exists():
        print("The MGF folder does not exist!")
        return
    if not glabel_path.exists():
        print("The gLabel folder does not exist!")
        return

    # Load mgf and gLabel files, then generate msp files for the corresponding names.
    # Iterate all the files in the MGF folder
    for mgf_file in mgf_path.iterdir():
        # Judge whether mgf is a folder, only open it as a mgf file
        if not mgf_file.is_file():
            print("There is a folder in the MGF folder!")
            continue
        if mgf_file.suffix != '.mgf':
           print("There is a file whose type is not mgf!")
           continue
        time_read_mgf = time.time()
        # Get the filename without the extension, mgf_file is the full "path/name" for the mgf file.
        mgf_stem = mgf_file.stem
        mgf_name = mgf_file.name
        # Generate the name for a glabel file, and the whole path/file for it.
        # add '-glabel' at the end for mgf_file to generate gLabel files
        # Example: 202110_Palleon_Cyno_In_vitro-20211027_Palleon_A_HILIC_HCDFT.mgf-glabel.txt
        glabel_name = f'{mgf_name}-glabel.txt'
        glabel_path_file = glabel_path / glabel_name
        # Generate the name for a msp file, and the whole path/file for it.
        msp_name = f'{mgf_stem}.msp'
        msp_path_file = msp_path / msp_name

        num_samples_one_file = mgf_glabel_to_msp(mgf_file, glabel_path_file, msp_path_file)

        total_num_files += 1
        total_num_samples += num_samples_one_file
        time_write_msp = time.time()
        process_time = time_write_msp - time_read_mgf
        print(f"The time for reading the mgf and glabel file, then write the msp file: {process_time} Seconds")

    time_end_mgf =  time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60
    total_time_hours = total_time_minutes / 60
    print(f"The total number of files is: {total_num_files}")
    print(f"The total number of samples is: {total_num_samples}")
    print(f"The end time is: {time_end_mgf}")
    print(f"The total time is: {total_time_seconds} Seconds, {total_time_minutes} Minutes, {total_time_hours} Hours.")


#  CLI (command line interface) for the input and output
def user_interface(input_type, input_path):
    """
    Get the folder by the input type and input path for the MSP files
    :param input_type: A string for the folder of input type;
    :param input_path: A string for the folder of input path.
    :return: The output of MSP files.
    """
    #  python construct_annotated_library.py -IT=N -IP=/data/Training-01-Human-285
    parse_mgf_files(input_type, input_path)

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


if __name__ == '__main__':
    # Please run the script with the following input format in Linux/Unix/Mac such as:
    # python construct_annotated_library.py -IT=N -IP=/data/Training-01-Human-285
    # Or Windows such as:
    # python construct_annotated_library.py -IT=O -IP=D:\\data\\Training-01-Human-285
    # For huge files with several hours, should use "nohup" and "&" to run in the background. For example:
    # nohup python -u construct_annotated_library.py -IT=N -IP=/data/Training-2021-D-Va--Human > D-Va–Human_out.out 2>&1 &

    try:
        user_interface(args.InputType, args.InputPath)
    except Exception as e:
        print(e)

