#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
########################################################################################################################
This script read csv file with glycopeptide sequence, then run the deep learning model for the predictions.
The output is a msp file similar to the glycopeptide spectral library.

Created on 01 June 2022.
Modified on 27 Aug 2022, to support 9 different types ions with maximum charge of 6.
########################################################################################################################
"""
__author__ = 'ZLiang'

from models import TestModel
from trainer import Trainer
import os
import torch as t
import argparse
from pathlib import Path
import time
import numpy as np
import pandas as pd
import spectral_library.parse_msp_to_csv as msp_to_csv
import spectral_library.calculate_fragment_mz as calculate_mz

os.environ['CUDA_VISIBLE_DEVICES'] = '1'
print(f"CUDA_VISIBLE_DEVICES: {os.environ['CUDA_VISIBLE_DEVICES']}")


class Config:
    # learning rate might be adjusted with the square root of (new_batch_size / old_batch_size).
    # Or in practical, adjusted with the ratio of (new_batch_size / old_batch_size).
    # Here for the first round, we use the ratio for the adjustment.
    lr = 0.0001
    batch_size = 8
    #batch_size = 16
    #batch_size = 32
    #batch_size = 64
    # batch_size = 128
    max_epoch = 500
    #gpu = False
    if t.cuda.is_available():
        gpu = True
    else:
        gpu = False

    print(f"The learning rate is: {lr}")
    print(f"The batch size is: {batch_size}")

opt = Config()

# Assume the maximum length of peptide is 32, and maximum length of glycan is 18.
# The maximum charge of glycopeptide is 6.
# The maximum number of ions is 3 (b, y, Y).
MAX_PEPTIDE_LENGTH = 32
MAX_GLYCAN_LENGTH = 18
MAX_NUM_CHARGES = 4
MAX_NUM_IONS = 9

# sequence length
SEQ_LEN = MAX_PEPTIDE_LENGTH + MAX_GLYCAN_LENGTH

# Change input_dim=24 to 26
net = TestModel(input_dim=26,
                n_tasks=MAX_NUM_IONS * MAX_NUM_CHARGES,
                embedding_dim=256,
                hidden_dim_lstm=128,
                hidden_dim_attention=32,
                n_lstm_layers=2,
                n_attention_heads=8,
                gpu=opt.gpu)
trainer = Trainer(net, opt)


# Predict the spectra from the input path, based on the saved models in saved_model.
def predict_by_sequence(input_path, saved_model, output_path):
    """
    :param input_path: the csv files for the training;
    :param saved_model: the pth file of the saved model for the prediction;
    :param output_path: saved msp files after prediction.
    :return:
    """
    input_folder = Path(input_path)
    # Check whether path name is a folder
    assert input_folder.is_dir(), "Input training path is wrong!"
    saved_pth = Path(saved_model)
    assert saved_pth.is_file(), "Input saved model is wrong!"
    output_folder = Path(output_path)
    assert output_folder.is_dir(), "Output saved model path is wrong!"
    
    trainer.load(saved_pth)
    print(f"The saved model is: {saved_pth}")

    time_start_predict = time.asctime(time.localtime(time.time()))
    start_time = time.time()
    print(f"The start time for predicting files: {time_start_predict}")

    num_train_files = 0
    num_total_samples = 0
    # Iterate all the pickle files in the input folder
    for input_file in input_folder.iterdir():
        print(f"The file for the prediction: {input_file}")

        num_train_files += 1
        df_sequence = pd.read_csv(input_file, sep='\t', header=0, encoding='utf8')
        num_total_samples += len(df_sequence)

        start_predict_time = time.time()

        output_dl_file = input_file.stem + "_dl" + str(".msp")
        output_dl_path = output_folder / output_dl_file
        output_db_file = input_file.stem + "_db" + str(".msp")
        output_db_path = output_folder / output_db_file

        with open(output_dl_path, 'w') as msp_dl_writer, open(output_db_path, 'w') as msp_db_writer:

            Xs = []
            X_metas = []

            for index, row in df_sequence.iterrows():
                peptide_name = df_sequence.iloc[index]['peptide']

                # Use 'X' to represent 'I' and 'L' that with the same masses,
                # save room for the fifth monosaccharide.
                peptide_name = peptide_name.replace('I', 'X')
                peptide_name = peptide_name.replace('L', 'X')
                # Assume the maximum length of peptide is 32, pad the left locations with "Z"
                left_peptide_length = MAX_PEPTIDE_LENGTH - len(peptide_name)
                peptide_pad_name = peptide_name + 'Z' * left_peptide_length

                monosaccharides = df_sequence.iloc[index]['glycan_denovo']
                # Assume the maximum length of glycan is max_glycan_length, pad the left locations with "Z"
                left_glycan_length = MAX_GLYCAN_LENGTH - len(monosaccharides)

                for monosaccharide, monosaccharide_code in msp_to_csv.monosaccharide_component_replacements.items():
                    # Replace monosaccharide with our special codes
                    monosaccharides = monosaccharides.replace(monosaccharide, monosaccharide_code)
                glycan_pad_name = monosaccharides + 'Z' * left_glycan_length
                glycopeptide_name = peptide_pad_name + glycan_pad_name
                glycopeptide_one_hot_encoded = msp_to_csv.one_hot_encode(glycopeptide_name,
                                                  msp_to_csv.amino_acid_monosaccharide_zero_codes)
                ion_charge = int(df_sequence.iloc[index]['charge'])

                Xs.append(np.array(glycopeptide_one_hot_encoded))
                X_metas.append(ion_charge)

            # Construct data frame for the deep learning model
            sample_input = [(X, X_meta, np.zeros((SEQ_LEN - 1, MAX_NUM_IONS * MAX_NUM_CHARGES))) for X, X_meta in zip(Xs, X_metas)]
            sample_pred = trainer.predict(sample_input)

            for i, pred in enumerate(sample_pred):
                glycopeptide_seq = [list(X) for X in sample_input[i][0]]
                charge = sample_input[i][1]
                glycopeptide_name = calculate_mz.CalculateFragmentMZ.reverse_one_hot_encode(
                    glycopeptide_seq, msp_to_csv.amino_acid_monosaccharide_zero_codes)

                # Get the peptide length and glycan length
                peptide_len = 0
                for c in glycopeptide_name[0: MAX_PEPTIDE_LENGTH - 1]:
                    peptide_len += 1
                    if c == 'Z':
                        break

                glycan_len = 0
                for c in glycopeptide_name[MAX_PEPTIDE_LENGTH: MAX_PEPTIDE_LENGTH + MAX_GLYCAN_LENGTH - 1]:
                    glycan_len += 1
                    if c == 'Z':
                        break

                #print("peptide_len=", peptide_len)
                #print("glycan_len=", glycan_len)

                calculator_mz = calculate_mz.CalculateMZ()
                for key, value in calculate_mz.dumb_reversal.items():
                    glycopeptide_name = glycopeptide_name.replace(key, value)

                mz = round(calculator_mz.calculate_mass(glycopeptide_name, charge), 4)
                mzs = []
                chemical_total_formulas=[]
                intensities = []
                ion_types = []
                ion_numbers = []
                positions = []
                ion_charges = []

                # the threshold for the output
                peaks = np.where(pred > 0.001)
                #print(peaks)
                # Consider four different types of ions: b ion, y ion, Y ion.
                # The max position of the ion is limited to "seq_len - 1".
                # For the length glycopeptide of 50, we have 49 for y.
                # Use "b, y, Y" to represent ions, and times the maximum charges, then the space is 9*4 = 36
                for position, peak_type in zip(*peaks):
                    # Convert peak_type back to ion_charge and ion_type
                    ion_charge = int(peak_type // MAX_NUM_IONS + 1)
                    ion_number = int(peak_type % MAX_NUM_IONS)
                    # ion_type should be 0, 1, 2, ... , 8 which represent b, b$, b-N(1), y, y$, y-N(1), Y, Y0, Y$ ions.
                    # position starts from 0, so should minus 1, and fragmentation should minus 1.
                    # Should consider the possibilities of different ions appearing to the positions
                    if ion_number == 0 and position < peptide_len - 1 - 1:
                        # b ion
                        ion = 'b'
                        pos = int(position + 1)
                    elif ion_number == 1 and position < peptide_len - 1 - 1:
                        # b$
                        ion = 'b$'
                        pos = int(position + 1)
                    elif ion_number == 2 and position < peptide_len - 1 - 1:
                        # b-N(1)
                        ion = 'b-N(1)'
                        pos = int(position + 1)
                    elif ion_number == 3 and position < peptide_len - 1 - 1:
                        # y ion
                        ion = 'y'
                        pos = int(MAX_PEPTIDE_LENGTH - position - 1)
                    elif ion_number == 4 and position < peptide_len - 1 - 1:
                        # y$
                        ion = 'y$'
                        pos = int(MAX_PEPTIDE_LENGTH - position - 1)
                    elif ion_number == 5 and position < peptide_len - 1 - 1:
                        # y-N(1)
                        ion = 'y-N(1)'
                        pos = int(MAX_PEPTIDE_LENGTH - position - 1)
                    elif ion_number == 6 and position == MAX_PEPTIDE_LENGTH - 1:
                        # Y0
                        ion = 'Y0'
                        pos = int(position + 1)
                    elif ion_number == 7 and position == MAX_PEPTIDE_LENGTH - 1:
                        # Y$
                        ion = 'Y$'
                        pos = int(position + 1)
                    elif ion_number == 8 and MAX_PEPTIDE_LENGTH - 1 < position < MAX_PEPTIDE_LENGTH + glycan_len - 1:
                        # Y ion
                        ion = 'Y'
                        pos = int(position + 1)
                    else:
                        continue
                    ion_types.append(ion)
                    ion_numbers.append(ion_number)
                    positions.append(pos)
                    ion_charges.append(ion_charge)

                    calculator_frag_mz = calculate_mz.CalculateFragmentMZ()
                    mz, chemical_formulas = calculator_frag_mz.get_frag_mz(glycopeptide_seq, pos, ion_number, ion_charge)

                    mzs.append(round(mz, 4))
                    chemical_total_formulas.append(chemical_formulas)
                    intensities.append(pred[position, peak_type])

                peptide = df_sequence.iloc[i]['peptide']
                glycan_denovo = df_sequence.iloc[i]['glycan_denovo']
                charge = str(df_sequence.iloc[i]['charge'])

                msp_dl_writer.write('Name: ' + peptide + '\t' + glycan_denovo + '\t' + charge + '\n')
                msp_dl_writer.write(f'MW: {str(mz)}' + '\n')
                msp_dl_writer.write('Comment: ' + 'CHARGE=' + str(charge) + '\t' +
                                 'GLYCOPEPTIDE=' + str(glycopeptide_name) + '\n')
                msp_dl_writer.write(f'Num Peaks: {len(ion_types)}' + '\n')

                msp_db_writer.write('Name: ' + peptide + '\t' + glycan_denovo + '\t' + charge + '+' + '\n')
                msp_db_writer.write(f'MW: {str(mz)}' + '\n')
                glycan_mono = "HNFAG"
                glycan_number = [0, 0, 0, 0, 0]
                for c in glycan_denovo:
                    glycan_number[glycan_mono.find(c)] += 1
                msp_db_writer.write('Comment: ' + 'CHARGE=' + str(charge) + '+' + '\t' + 'PEPTIDE=' + peptide + '\t'
                                    + 'GLYCAN(H,N,F,A,G)=' + str(glycan_number)[1:-1].replace(" ", "") + '\n')
                msp_db_writer.write(f'Num Peaks: {len(ion_types)}' + '\n')

                for j in range(len(ion_types)):
                    msp_dl_writer.write(str(mzs[j]) + '\t' + str(intensities[j]) + '\t' + str(ion_types[j]) + str(positions[j])
                                        + '+' + str(ion_charges[j]) + '\t' + chemical_total_formulas[j] + '\n')

                    ion_annotation = calculator_frag_mz.get_frag_annotation(glycopeptide_name,
                                                                            positions[j], ion_numbers[j], ion_charges[j])
                    msp_db_writer.write(str(mzs[j]) + '\t' + str(intensities[j]) + '\t' + ion_annotation +
                                        '\t' + chemical_total_formulas[j] +'\n')
                msp_dl_writer.write('\n')
                msp_db_writer.write('\n')

            end_predict_time = time.time()
            predict_time_seconds = round(end_predict_time - start_predict_time, 2)
            predict_time_minutes = round(predict_time_seconds / 60, 2)
            predict_time_hours = round(predict_time_minutes / 60, 2)
            print(f"Time for prediction: {predict_time_seconds} S, {predict_time_minutes} M, {predict_time_hours} H.")
            print(f"Successfully generate a prediction deep learning msp file: {output_dl_file}")
            print(f"Successfully generate a prediction database msp file: {output_db_file}")

    print("The total number of training files: %d" % num_train_files)
    print("The total number of samples in the training datasets: %d" % num_total_samples)

    time_end_predict = time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = round(end_time - start_time, 2)
    total_time_minutes = round(total_time_seconds / 60, 2)
    total_time_hours = round(total_time_minutes / 60, 2)

    print(f"The end time is: {time_end_predict}")
    print(f"The total time is: {total_time_seconds} Seconds, {total_time_minutes} Minutes, {total_time_hours} Hours.")


#  CLI (command line interface) for the input and output
def user_interface(input_path, saved_model, output_path):
    """
    :param input_path: path for the input csv files;
    :param saved_model: the saved model for the prediction;
    :param output_path: path for the output msp files;
    :return:
    """
    # python predict_by_sequence.py -IP=/data/input_sequence
    #   -SM=/data/saved_models/train-345-76-Full-top-1_epoch-50_batch-64_lr-4_24-May-2022/model-49.pth
    #   -OP=D:\data\predict_spectra
    predict_by_sequence(input_path, saved_model, output_path)


"""
Input parameters for the user interface:
    1   Input Training Path    
    2   Input Saved Model 
    2   Output Saved Model Path
"""
parser = argparse.ArgumentParser(description='Input parameters to run the script.')
parser.add_argument('--InputType', '-IT',
                    help='Input File parameter，required, no default. such as /data/testing',
                    required=False)
parser.add_argument('--TopNumber', '-TN',
                    help='Input Saved Model parameter，required，no default. Such as '
                         '/data/saved_models/train-345-76-Full-top-1_epoch-50_batch-32_lr-2_26-May-2022/model-49.pth.',
                    required=False)
parser.add_argument('--InputPath', '-IP',
                    help='Output Path parameter，required，no default. Such as '
                         '/data/prediction/',
                    required=False)

args = parser.parse_args()


if __name__ == "__main__":
    # Please run the script with the following input format in Linux/Unix/Mac such as:
    # python predict_by_sequence.py -IP=/data/input_sequence
    #  -SM=/data/saved_models/train-345-76-Full-top-1_epoch-50_batch-32_lr-2_26-May-2022/model-49.pth
    #  -OP=/data/predict_spectra
    #
    # python predict_by_sequence.py -IP=D:\data\input_sequence
    #  -SM=D:\data\saved_models\train-345-76-Full-top-1_epoch-50_batch-32_lr-2_26-May-2022\model-49.pth
    #  -OP=D:\data\predict_spectra
    #
    # For huge files with several hours, should use "nohup" and "&" to run in the background. For example:
    # nohup python -u predict_by_sequence.py -IP=/data/input_sequence
    #  -SM=/data/saved_models/train-345-76-Full-top-1_epoch-50_batch-32_lr-2_26-May-2022/model-49.pth
    #  -OP=/data/predict_spectra
    #  > Training-01-Human-285_Testing-01_Batch-8.out 2>&1 &

    try:
        user_interface(args.InputType, args.TopNumber, args.InputPath)
    except Exception as e:
        print(e)