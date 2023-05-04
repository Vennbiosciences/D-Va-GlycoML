#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
################################################################################
This script combines mgf, gGlyco3, and gLabel to generate spectral libraries.
This spectral library only keep the annotated b-ions, y-ions, Y-ions, etc, and
ignore the unannotated ones.

For the purpose of machine learning and deep learning.
Created on 10 August 2021, modified on 17 September 2021.

TBD: Need to modified for keeping all the information.
################################################################################
"""
__author__ = 'ZLiang'


from pyteomics import mgf, auxiliary
import pandas as pd
import numpy as np
import decimal
import os

# define binary search for the glabel['spec']
def binary_search(spectra, target):
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

# define file path, and mgf, msp, pglyco, glabel files.
# Should consider all the folds in file_path, to generate the library automatically. And there is a raw file folder,
# which we need to ignore it.

# file_path = '/home/zliang/orbitrap_data/data/21-08-IO/'
msp_path = '/home/zliang/orbitrap_data/Spectral_Libraries/Annotated_Spectrum_Prediction/DDA_2021/MSP/'
fold_path = '/home/zliang/orbitrap_data/data_fragmentation_training/'

file_path = fold_path + '210714-50cm-column/'
# file_path = fold_path + '210723-FAIMS-Test/'
# file_path = fold_path + '210714-50cm-column/'
# file_path = fold_path + '210714-50cm-column/'
# file_path = fold_path + '210714-50cm-column/'

mgf_path = file_path + 'MGF/'
glabel_path = file_path + 'gLabel-ions/'

# Load mgf and pLabel files, generate msp files
# Explore all the files in MGF folder
if os.path.exists(mgf_path):
    print('Good')
    # get all the mgf files
    mgf_files = os.listdir(mgf_path)
    for mgf_file in mgf_files:
        # Judge whether mgf is a folder, only open it as a mgf file
        if not os.path.isdir(mgf_file):
            # get the filename without the extension using 'rsplit',
            mgf_name = mgf_file.rsplit('.', 1)[0]
            # remove 'formatted' by '_'
            file_name = mgf_name.split('_')[0]
            # add 'N-' at the start for pGlyco
            pGlyco_name = 'N-' + file_name
            # add '.msp' to generate the msp file
            glabel_file = glabel_path + mgf_file + '-glabel.txt'
            # Load gLabel
            glabel = pd.read_csv(glabel_file, sep='\t', header=0, encoding='utf8', engine='python',
                                 error_bad_lines=False)
            glabel_sort = glabel.sort_values(by=["spec"])
            # create msp file
            msp_name = pGlyco_name + '.msp'
            msp_file = msp_path + msp_name
            mgf_path_file = mgf_path + mgf_file
            with open(msp_file, 'w') as writer:
                with mgf.read(mgf_path_file) as reader:
                    for spectrum in reader:
                        index = binary_search(glabel_sort['spec'], spectrum['params']['title'])
                        if index != -1:
                            anno_array = ['' for i in range(len(spectrum['m/z array']))]
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
                            writer.write('Name: ' + spectrum['params']['title'] + '\n')
                            writer.write('MW: ' + str(spectrum['params']['pepmass'][0]) + '\n')
                            writer.write('Comment: ' + 'CHARGE=' + str(spectrum['params']['charge']) + '\t' + \
                                         'RTINSECONDS=' + str(spectrum['params']['rtinseconds']) + '\t' + \
                                         'PEPTIDE=' + glabel_sort.iloc[index]['peptide'] + '\t' + \
                                         'MODIFICATIONS=' + str(glabel_sort.iloc[index]['modinfo']) + '\t' + \
                                         'GLYCAN(H,N,F,A)=' + str(glabel_sort.iloc[index]['glycan(H,N,F,A)']) + '\t' + \
                                         'GLYCANS=' + str(glabel_sort.iloc[index]['formula']) + '\n')
                            writer.write('Num Peaks: ' + str(len(spectrum_mz)) + '\n')
                            for i in range(len(spectrum_mz)):
                                writer.write(
                                    str(spectrum_mz[i]) + '\t' + str(spectrum_intensity[i]) + '\t' \
                                    + str(spectrum_annotation[i]) + '\n')
                            writer.write('\n')
                writer.close()
        else:
            print('The mgf file you specified is a folder. \n')

