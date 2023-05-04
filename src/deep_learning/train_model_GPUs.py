#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
################################################################################
This script train the deep learning model for the predictions.

Created on 5 April 2022.
Modified on 19 April 2022, for the command line interface and batch processing.
################################################################################
"""
__author__ = 'ZLiang'

import numpy as np
import pandas as pd
from biLSTM import BiLSTM, MultiheadAttention
from models import TestModel
from trainer_GPUs import Trainer_GPUs
from sklearn.metrics import precision_score, recall_score
import os
import pickle
import argparse
from pathlib import Path
import torch as t
import time

os.environ['CUDA_VISIBLE_DEVICES'] = '3'
#os.environ['CUDA_VISIBLE_DEVICES'] = '0, 1, 2, 3'
print(os.environ['CUDA_VISIBLE_DEVICES'])


class Config:
    # learning rate might be adjusted with the square root of (new_batch_size / old_batch_size).
    # Or in practical, adjusted with the ratio of (new_batch_size / old_batch_size).
    lr = 0.0001
    batch_size = 8
    #batch_size = 16
    #batch_size = 32
    #batch_size = 64
    # batch_size = 128
    max_epoch = 50
    # gpu = False
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
trainer_gpus = Trainer_GPUs(net, opt)


# Training the model from the train_path, testing from the test_path, save the trained models in saved_model_path.
def train_model(train_path, test_path, saved_model_path):
    """
    :param train_path: the pkl files for the training;
    :param test_path: the pkl files for the testing;
    :param saved_model_path: saved models after training.
    :return:
    """
    #model_saved = "/data/saved_models/model_07_Mar.pth"
    # trainer.load(model_saved)

    train_folder = Path(train_path)
    # Check whether path name is a folder
    assert train_folder.is_dir(), "Input training path is wrong!"
    test_folder = Path(test_path)
    assert test_folder.is_dir(), "Input testing path is wrong!"
    saved_model_folder = Path(saved_model_path)
    assert saved_model_folder.is_dir(), "Output saved model path is wrong!"

    time_start_train = time.asctime(time.localtime(time.time()))
    start_time = time.time()
    print(f"The start time for training files: {time_start_train}")

    train_input = []

    num_train_files = 0
    # Iterate all the pickle files in the training folder
    for train_file in train_folder.iterdir():
        num_train_files += 1
        train_pickle = pickle.load(open(train_file, 'rb'))
        # add each pickle file into the list, should use "extend", not "append".
        train_input.extend(train_pickle)

    print("The total number of training files: %d" % num_train_files)
    print("The total number of samples in the training datasets: %d" % len(train_input))

    """
    # Get the total number of files the in test folder
    num_test_files = sum(1 for _ in test_folder.iterdir())
    print("The total number of testing files: %d" % num_test_files)

    test_file_lists = [[] for _ in range(num_test_files)]
    test_name_lists: List[List[str]] = [[] for _ in range(num_test_files)]
    # Iterate all the pickle files in the testing folder
    i: int
    for i, test_file in enumerate(test_folder.iterdir()):
        print(i)
        test_pickle = pickle.load(open(test_file, 'rb'))
        test_file_lists[i] = test_pickle
        test_name_lists[i] = test_file.stem

    cos_sims_last_lists = [[0] for _ in range(num_test_files)]
    cos_sims_current_lists = [[0] for _ in range(num_test_files)]
    sims_difference_lists = [[0] for _ in range(num_test_files)]
    precisions_current_lists = [[0] for _ in range(num_test_files)]
    recalls_current_lists = [[0] for _ in range(num_test_files)]

    for i in range(100):
        start_fold_time = time.time()
        print("Start fold %d" % i)

        # Calculate the Cos Similarity difference of two continuous folds
        if i > 0:
            for j in range(len(test_name_lists)):
                print(f"The number of samples in the testing file of {test_name_lists[j]}: {len(test_file_lists[j])}")
                precisions_current_lists[j], recalls_current_lists[j], cos_sims_current_lists[j] = \
                    evaluate(test_file_lists[j], trainer.predict(test_file_lists[j]))

        if i > 1:
            for j in range(len(test_name_lists)):
                sims_difference_lists[j] = abs(np.mean(cos_sims_current_lists[j]) - np.mean(cos_sims_last_lists[j]))
                print(f"The difference of Cosine Similarity for {test_name_lists[j]}: {sims_difference_lists[j]}")

        for j in range(len(test_name_lists)):
            cos_sims_last_lists[j] = cos_sims_current_lists[j]

        end_fold_time = time.time()
        fold_time = end_fold_time - start_fold_time

        print("Time for evaluation:", fold_time)

        """
    for i in range(40):
        print("Start fold %d" % i)

        start_train_time = time.time()
        trainer_gpus.train(train_input, n_epochs=5)
        end_train_time = time.time()

        train_time_seconds = end_train_time - start_train_time
        train_time_minutes = train_time_seconds / 60
        train_time_hours = train_time_minutes / 60

        print(f"Time for training: {train_time_seconds} S, {train_time_minutes} M, {train_time_hours} H.")

        trainer_gpus.save(os.path.join(saved_model_folder, 'model-%d.pth' % i))

        end_save_time = time.time()
        save_time = end_save_time - end_train_time

        print(f"Time for saving: {save_time} seconds")

    time_end_train = time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60
    total_time_hours = total_time_minutes / 60

    print(f"The end time is: {time_end_train}")
    print(f"The total time is: {total_time_seconds} Seconds, {total_time_minutes} Minutes, {total_time_hours} Hours.")


#  CLI (command line interface) for the input and output
def user_interface(train_path, test_path, saved_model_path):
    """
    :param train_path: path for the training
    :param test_path: path for the testing
    :param saved_model_path: path for the saved model
    :return:
    """
    #  python train_model_GPUs.py -TA=/data/Training-01-Human-285/N-GP-PKL
    #       -TE=/data/Testing-01-Different-HCD-energies-24/Energy-01-HCD-15-20-34-37-40/N-GP-PKL
    #       -SM=/data/saved_models/train-15_test-4_batch-8_14-Apr-2022
    train_model(train_path, test_path, saved_model_path)


"""
Input parameters for the user interface:
    1   Input Training Path
    2   Input Testing Path
    3   Output Saved Model Path
"""
parser = argparse.ArgumentParser(description='Input parameters to run the script.')
parser.add_argument('--TrainingPath', '-TA',
                    help='Input Training Path parameter，required, no default. such as '
                         '/data/Training-01-Human-285/N-GP-PKL.',
                    required=False)
parser.add_argument('--TestingPath', '-TE',
                    help='Input Testing Path parameter，required，no default. Such as '
                         '/data/Testing-01-Different-HCD-energies-24/Energy-01-HCD-15-20-34-37-40/N-GP-PKL.',
                    required=False)
parser.add_argument('--SavedModelPath', '-SM',
                    help='Output Saved Model Path parameter，required，no default. Such as '
                         '/data/saved_models/train-15_test-4_batch-8_14-Apr-2022.',
                    required=False)

args = parser.parse_args()

if __name__ == "__main__":
    # Please run the script with the following input format in Linux/Unix/Mac such as:
    # python train_model_GPUs.py -TA=/data/Training-01-Human-285/N-GP-PKL
    # -TE=/data/Testing-01-Different-HCD-energies-24/Energy-01-HCD-15-20-34-37-40/N-GP-PKL
    # -SM=/data/saved_models/train-15_test-4_batch-8_14-Apr-2022
    #
    # python train_model_GPUs.py -TA=D:\data\Training-01-Human-285\N-GP-PKL
    # -TE=D:\data\Testing-01-Different-HCD-energies-24\Energy-01-HCD-15-20-34-37-40\N-GP-PKL
    # -SM=D:\data\saved_models\train-15_test-4_batch-8_14-Apr-2022
    #
    # For huge files with several hours, should use "nohup" and "&" to run in the background. For example:
    # nohup python -u train_model_GPUs.py -TA=/data/Training-01-Human-285/N-GP-PKL
    # -TE=/data/Testing-01-Different-HCD-energies-24/Energy-01-HCD-15-20-34-37-40/N-GP-PKL
    # -SM=/data/saved_models/train-15_test-4_batch-8_14-Apr-2022
    # > Training-01-Human-285_Testing-01_Batch-8.out 2>&1 &

    try:
        user_interface(args.TrainingPath, args.TestingPath, args.SavedModelPath)
    except Exception as e:
        print(e)