#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
########################################################################################################################
This script train the deep learning model for the predictions.
By reading all the files in training folder into the memory, this approach will train all the datasets for one epoch.

Created on 23 March 2022.
Modified on 14 April 2022, for the command line interface and batch processing.
Modified on 21 April 2022, lock down for the batch process, by reading all the training files into memory.
Modified on 22 April 2022, focus on the training dataset of Training-01-Human-285 for top one de novo candidate.
Modified on 18 May 2022, output information for each epoch, and the maximum number of epochs is 50, for cyno and mouse.
########################################################################################################################
"""
__author__ = 'ZLiang'

from models import TestModel
from trainer import Trainer
import os
import pickle
import torch as t
import argparse
from pathlib import Path
import time

os.environ['CUDA_VISIBLE_DEVICES'] = '0'
print(f"CUDA_VISIBLE_DEVICES: {os.environ['CUDA_VISIBLE_DEVICES']}.")


class Config:
    # learning rate might be adjusted with the square root of (new_batch_size / old_batch_size).
    # Or in practical, adjusted with the ratio of (new_batch_size / old_batch_size).
    # Here for the first round, we use the ratio for the adjustment.
    lr = 0.0008
    #batch_size = 8
    #batch_size = 16
    #batch_size = 32
    batch_size = 64
    #batch_size = 128
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


# Training the model from the train_path, save the trained models in saved_model_path.
def train_model_cyno(train_path, saved_model_path):
    """
    :param train_path: the pkl files for the training;
    :param saved_model_path: saved models after training.
    :return:
    """
    train_folder = Path(train_path)
    # Check whether path name is a folder
    assert train_folder.is_dir(), "Input training path is wrong!"
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

    for i in range(50):
        print("Start epoch %d" % i)
        start_train_time = time.time()

        trainer.train(train_input, n_epochs=1)

        end_train_time = time.time()
        train_time_seconds = round(end_train_time - start_train_time, 2)
        train_time_minutes = round(train_time_seconds / 60, 2)
        train_time_hours = round(train_time_minutes / 60, 2)
        print(f"Time for training: {train_time_seconds} S, {train_time_minutes} M, {train_time_hours} H.")

        trainer.save(os.path.join(saved_model_folder, 'model-%d.pth' % i))

        end_save_time = time.time()
        save_time = round(end_save_time - end_train_time, 2)

        print(f"Time for saving: {save_time} seconds")

    time_end_train = time.asctime(time.localtime(time.time()))
    end_time = time.time()
    total_time_seconds = round(end_time - start_time, 2)
    total_time_minutes = round(total_time_seconds / 60, 2)
    total_time_hours = round(total_time_minutes / 60, 2)

    print(f"The end time is: {time_end_train}")
    print(f"The total time is: {total_time_seconds} Seconds, {total_time_minutes} Minutes, {total_time_hours} Hours.")


#  CLI (command line interface) for the input and output
def user_interface(train_path, saved_model_path):
    """
    :param train_path: path for the training
    :param saved_model_path: path for the saved model
    :return:
    """
    #  python train_model_cyno.py -TA=/data/Training-01-Human-285/N-GP-PKL
    #       -SM=/data/saved_models/train-15_test-4_batch-8_14-Apr-2022
    train_model_cyno(train_path, saved_model_path)


"""
Input parameters for the user interface:
    1   Input Training Path    
    2   Output Saved Model Path
"""
parser = argparse.ArgumentParser(description='Input parameters to run the script.')
parser.add_argument('--TrainingPath', '-TA',
                    help='Input Training Path parameter，required, no default. such as '
                         '/data/Training-01-Human-285/N-GP-PKL.',
                    required=False)
parser.add_argument('--SavedModelPath', '-SM',
                    help='Output Saved Model Path parameter，required，no default. Such as '
                         '/data/saved_models/train-15_test-4_batch-8_14-Apr-2022.',
                    required=False)

args = parser.parse_args()


if __name__ == "__main__":
    # Please run the script with the following input format in Linux/Unix/Mac such as:
    # python train_model_cyno.py -TA=/data/Training-01-Human-285/N-GP-PKL
    #  -SM=/data/saved_models/train-15_test-4_batch-8_14-Apr-2022
    #
    # python train_model_cyno.py -TA=D:\data\Training-01-Human-285\N-GP-PKL
    #  -SM=D:\data\saved_models\train-15_test-4_batch-8_14-Apr-2022
    #
    # For huge files with several hours, should use "nohup" and "&" to run in the background. For example:
    # nohup python -u train_model_cyno.py -TA=/data/Training-01-Human-285/N-GP-PKL
    # -SM=/data/saved_models/train-15_test-4_batch-8_14-Apr-2022
    # > Training-01-Human-285_Testing-01_Batch-8.out 2>&1 &

    try:
        user_interface(args.TrainingPath, args.SavedModelPath)
    except Exception as e:
        print(e)