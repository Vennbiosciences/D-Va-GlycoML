#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
########################################################################################################################
This script train the deep learning model for the predictions.

Created on 23 March 2022.
Modified on 14 April 2022, for the command line interface and batch processing.
Modified on 27 April 2022, for the batch processing of testing and outputting the prediction for each spectrum.
########################################################################################################################
"""
__author__ = 'ZLiang'

from typing import List, Any

import numpy as np
import pandas as pd
from biLSTM import BiLSTM, MultiheadAttention
from models import TestModel
from trainer import Trainer
from sklearn.metrics import precision_score, recall_score
import os
import pickle
import torch as t
import argparse
from pathlib import Path
import time

os.environ['CUDA_VISIBLE_DEVICES'] = '1'
#os.environ['CUDA_VISIBLE_DEVICES'] = '0, 1, 2, 3'
print(f"CUDA_VISIBLE_DEVICES: {os.environ['CUDA_VISIBLE_DEVICES']}")


class Config:
    lr = 0.0001
    #batch_size = 8
    #batch_size = 16
    batch_size = 32
    #batch_size = 64
    # batch_size = 128
    max_epoch = 50
    #gpu = False
    if t.cuda.is_available():
        gpu = True
    else:
        gpu = False


opt = Config()

# Assume the maximum length of peptide is 32, and maximum length of glycan is 18.
# The maximum charge of glycopeptide is 6.
# The maximum number of ions is 3 (b, y, Y).
MAX_PEPTIDE_LENGTH = 32
MAX_GLYCAN_LENGTH = 18
MAX_NUM_CHARGES = 6
MAX_NUM_IONS = 3


def evaluate(inputs, preds, thr=0.01):
    precisions = []
    recalls = []
    cos_sims = []

    for inp, pred in zip(inputs, preds):
        label = inp[2]
        _label = label > thr
        _pred = pred > thr
        precisions.append(precision_score(_label.flatten(), _pred.flatten()))
        recalls.append(recall_score(_label.flatten(), _pred.flatten()))
        sim = np.sum(label * pred) / np.sqrt(np.sum(label * label) * np.sum(pred * pred))
        cos_sims.append(sim)
    print("Precision for Mean %f" % np.mean(precisions))
    print("Precision for Median %f" % np.median(precisions))
    print("Recall for Mean %f" % np.mean(recalls))
    print("Recall for Median %f" % np.median(recalls))
    print("Cos Similarity for Mean %f" % np.mean(cos_sims))
    print("Cos Similarity for Median %f" % np.median(cos_sims))
    return precisions, recalls, cos_sims


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


# Training the model from the train_path, testing from the test_path, save the trained models in saved_model_path.
def train_model(model_file, input_type, top_number, test_path):
    """
    :param model_file: the pth file for the trained model;
    :param input_type: A string for the folder of input type, such as "N" or "O";
    :param top_number: A value for the top number of de novo sequencing, such as "10";
    :param test_path: A string for the input test folder, such as "/data/Training-01-Human-285".
    :return:
    """
    train_name = Path(model_file)
    # Check whether path name is a file
    assert train_name.is_file(), "Input model file is wrong!"
    path_name = Path(test_path)
    # Check whether test name is a folder
    assert path_name.is_dir(), "Input path is wrong!"

    top_number = int(top_number)
    assert top_number > 0 and top_number < 101, "Wrong Top Number!"

    if input_type == 'N':
        pkl_name = 'N-GP-PKL-TOP-' + str(top_number)
        pkl_path = path_name / pkl_name
        test_name = pkl_name + '-TEST'
        test_path = path_name / test_name
    elif input_type == 'O':
        pkl_name = 'O-GP-PKL-TOP-' + str(top_number)
        pkl_path = path_name / pkl_name
        test_name = pkl_name + '-TEST'
        test_path = path_name / test_name
    else:
        print("Error Input Type!")
        return

    time_start_test = time.asctime(time.localtime(time.time()))
    start_time = time.time()
    print(f"The start time for testing files: {time_start_test}")
    print(f"The pkl path is: {pkl_path.absolute()}")
    print(f"The test path is: {test_path.absolute()}")

    # Make dir for the test folder
    test_path.mkdir(parents=True, exist_ok=True)

    trainer.load(train_name)

    time_start_train = time.asctime(time.localtime(time.time()))

    print(f"The start time for training files: {time_start_train}")

    # Get the total number of files the in test folder
    num_test_files = sum(1 for _ in pkl_path.iterdir())
    print("The total number of testing files: %d" % num_test_files)

    test_file_lists = [[] for _ in range(num_test_files)]
    test_name_lists: List[List[str]] = [[] for _ in range(num_test_files)]
    # Iterate all the pickle files in the testing folder
    i: int
    for i, test_file in enumerate(pkl_path.iterdir()):
        test_pickle = pickle.load(open(test_file, 'rb'))
        test_file_lists[i] = test_pickle
        test_name_lists[i] = test_file.stem

    cos_sims_lists = [[0] for _ in range(num_test_files)]
    precisions_lists = [[0] for _ in range(num_test_files)]
    recalls_lists = [[0] for _ in range(num_test_files)]

    precisions_total = []
    recalls_total = []
    cos_sims_total = []

    # Predict the spectra with trained model for different testing files, and save the results to txt files.
    for j in range(len(test_name_lists)):
        print(f"The number of samples in the testing file of {test_name_lists[j]}: {len(test_file_lists[j])}")
        precisions_lists[j], recalls_lists[j], cos_sims_lists[j] = \
            evaluate(test_file_lists[j], trainer.predict(test_file_lists[j]))
        precisions_total.extend(precisions_lists[j])
        recalls_total.extend(recalls_lists[j])
        cos_sims_total.extend(cos_sims_lists[j])

        result_name = f"{test_name_lists[j]}.txt"
        result_file = test_path / result_name
        with open(result_file, 'w') as test_result_file:
            test_result_file.write(f"The total number of predicted samples in the testing file of "
                                   f"{test_name_lists[j]}: {len(test_file_lists[j])} \n")
            test_result_file.write(f"Precision for Mean {np.mean(precisions_lists[j])} \n")
            test_result_file.write(f"Precision for Median {np.median(precisions_lists[j])} \n")
            test_result_file.write(f"Recall for Mean {np.mean(recalls_lists[j])} \n")
            test_result_file.write(f"Recall for Median {np.median(recalls_lists[j])} \n")
            test_result_file.write(f"Cos Similarity for Mean {np.mean(cos_sims_lists[j])} \n")
            test_result_file.write(f"Cos Similarity for Median {np.median(cos_sims_lists[j])} \n")
            for k in range(len(cos_sims_lists[j])):
                test_result_file.write(f"{k+1} \t {cos_sims_lists[j][k]} \n")

    print("")
    print(f"The total number of predicted samples in all the testing file: {len(precisions_total)}")
    print(f"The Total Precision for Mean: {np.mean(precisions_total)}")
    print(f"The Total Precision for Median: {np.median(precisions_total)}")
    print(f"The Total Recall for Mean: {np.mean(recalls_total)}")
    print(f"The Total Recall for Median: {np.median(recalls_total)}")
    print(f"The Total Cos Similarity for Mean: {np.mean(cos_sims_total)}")
    print(f"The Total Cos Similarity for Median: {np.median(cos_sims_total)} \n")

    time_end_predict = time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60
    total_time_hours = total_time_minutes / 60
    print(f"The end time is: {time_end_predict}")
    print(f"The total time is: {total_time_seconds} Seconds, {total_time_minutes} Minutes, {total_time_hours} Hours.")


#  CLI (command line interface) for the input and output
def user_interface(model_file, input_type, top_number, test_path):
    """
    :param model_file: file for the trained model
    :param input_type: A string for the folder of input type;
    :param top_number: A value for the top number of de novo sequencing, such as "10";
    :param test_path: A string for the folder of input testing path.
    :return:
    """
    #  python predict_trained_model.py -MF=/data/saved_models/train-285-top-one_test-4_batch-32_lr-2_22-Apr-2022/train-285-top-one_test-4_batch-32_lr-2_22-Apr-2022_model-99.pth
    #   -IT=N -TN=1 -TP=/data/Testing-01-Different-HCD-energies-24/Energy-01-HCD-15-20-34-37-40
    train_model(model_file, input_type, top_number, test_path)

"""
Input parameters for the user interface:
    1   Input Model File
    2   Input Glyco Peptide Type
    3   Input Top Number of De Novo Sequencing
    4   Input and Output Testing Path
"""
parser = argparse.ArgumentParser(description='Input parameters to run the script.')
parser.add_argument('--ModelFile', '-MF',
                    help='Input Trained File parameter，required, no default. such as '
                         '/data/saved_models/train-285-top-one_test-4_batch-32_lr-2_22-Apr-2022',
                    required=False)
parser.add_argument('--InputType', '-IT',
                    help='Input Format parameter，required, has default. N-NlinkedGlycoPeptide; O-OlinkedGlycoPeptide.',
                    required=False, default='N')
parser.add_argument('--TopNumber', '-TN',
                    help='Top Number parameter，required, has default. The value of the top number for de novo '
                         'sequencing.', required=False, default='1')
parser.add_argument('--TestingPath', '-TP',
                    help='Input and Output Testing Path parameter，required，no default. Such as '
                         '/data/Testing-01-Different-HCD-energies-24/Energy-01-HCD-15-20-34-37-40/',
                    required=False)

args = parser.parse_args()


if __name__ == "__main__":
    # Please run the script with the following input format in Linux/Unix/Mac such as:
    # python predict_trained_model.py -MF=/data/saved_models/train-285-top-one_test-4_batch-32_lr-2_22-Apr-2022/train-285-top-one_test-4_batch-32_lr-2_22-Apr-2022_model-99.pth
    # -IT=N -TN=1 -TP=/data/Testing-01-Different-HCD-energies-24/Energy-01-HCD-15-20-34-37-40
    #
    # Or Windows such as:
    # python predict_trained_model.py -MF=d:\data\saved_models\train-285-top-one_test-4_batch-32_lr-2_22-Apr-2022\train-285-top-one_test-4_batch-32_lr-2_22-Apr-2022_model-99.pth
    # -IT=N -TN=1 -TP=D:\data\Testing-01-Different-HCD-energies-24\Energy-01-HCD-15-20-34-37-40
    #
    # For huge files with several hours, should use "nohup" and "&" to run in the background. For example:
    # nohup python -u predict_trained_model.py -MF=/data/saved_models/train-285-top-one_test-4_batch-32_lr-2_22-Apr-2022/train-285-top-one_test-4_batch-32_lr-2_22-Apr-2022_model-99.pth
    # -IT=N -TN=1 -TP=/data/Testing-01-Different-HCD-energies-24/Energy-01-HCD-15-20-34-37-40
    # > Training-01-Human-285_Testing-01_Batch-8.out 2>&1 &

    try:
        user_interface(args.ModelFile, args.InputType, args.TopNumber, args.TestingPath)
    except Exception as e:
        print(e)