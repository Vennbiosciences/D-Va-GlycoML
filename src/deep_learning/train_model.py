#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
################################################################################
This script train the deep learning model for the predictions.

Created on 26 October 2021.
Modified on 28 February 2022.
################################################################################
"""
__author__ = 'ZLiang'

import numpy as np
import pandas as pd
from biLSTM import BiLSTM, MultiheadAttention
from models import TestModel
from trainer import Trainer
from sklearn.metrics import precision_score, recall_score
import os
import pickle
from pathlib import Path
import torch as t
import time

os.environ['CUDA_VISIBLE_DEVICES'] = '0, 1, 2, 3'
print(os.environ['CUDA_VISIBLE_DEVICES'])

class Config:
    lr = 0.0001
    batch_size = 8
    #batch_size = 16
    #batch_size = 32
    #batch_size = 64
    #batch_size = 128
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
max_peptide_length = 32
max_glycan_length = 18
max_num_charges = 6
max_num_ions = 3

def evaluate(inputs, preds, thr=0.01):
    precisions = []
    recalls = []
    cos_sims = []

    for inp, pred in zip(inputs, preds):
        label = inp[2]
        #if opt.gpu:
        #     label = np.array(label)
        #     label = t.from_numpy(label)
        #     pred = np.array(pred)
        #    pred = t.from_numpy(pred)
        #    label = label.cuda()
        #    pred = pred.cuda()

        _label = label > thr
        _pred = pred > thr
        precisions.append(precision_score(_label.flatten(), _pred.flatten()))
        recalls.append(recall_score(_label.flatten(), _pred.flatten()))
        sim = np.sum(label * pred) / np.sqrt(np.sum(label * label) * np.sum(pred * pred))
        cos_sims.append(sim)
    print("Precision %f" % np.mean(precisions))
    print("Recall %f" % np.mean(recalls))
    print("Cos similarity %f" % np.mean(cos_sims))
    return precisions, recalls, cos_sims

"""
neutral_loss_choices = [0, 17, 18, 35, 36, 44, 46]
n_neutral_losses = len(neutral_loss_choices)
n_charges = 7


# df = pd.read_csv("../data/output.csv")
sptxt_file = "/home/zliang_venn_bio/spectranet_data/sptxt_output.csv"
csv_file = "/data/Training-2021-D-Va--Human/CSV/210908-RD-SS-21-09-Enrich-Lot2_8-02-1_0ul-check_HCDFT.csv"

#pickle_file ="/home/zliang/orbitrap_data/data/Glycan-DeNovo-Test/MSP/sptxt_output.pkl"

#df = pd.read_csv(csv_file, encoding = 'utf-8')
df = pd.read_csv(sptxt_file, encoding = 'utf-8')
# pickle_output = open(pickle_file,'wb')
# pickle.dump(df,pickle_output)
# pickle_output.close()


Xs = []
ys = []
X_metas = []
for line in np.array(df):
  '''
0  ['acetyl',
1 'name',
2 'charge',
3 'precursor',
4 'mz',
5 'intensity',
6 'ion',
7 'position',
8 'neutral_loss',
9 'ion_charge',
10 'delta']
  '''
  Xs.append(np.array(eval(line[1])))
  X_metas.append((line[0], line[2]))
  seq_len = np.array(eval(line[1])).shape[0]

  intensities = np.array(eval(line[5])).astype(float)
  ions = np.array(eval(line[6])).astype(int)
  neutral_losses = np.array(eval(line[8])).astype(int)
  ion_charges = np.array(eval(line[9])).astype(int)
  positions = np.array(eval(line[7])).astype(int)
  offsets = np.array(eval(line[10])).astype(float)
  assert np.min(positions) >= 1
  assert np.max(positions) <= seq_len
  y = np.zeros((seq_len - 1, 2*n_charges*n_neutral_losses))
  for i in range(len(intensities)):
    if np.abs(offsets[i]) > 0.1 or \
       ion_charges[i] > n_charges or \
       neutral_losses[i] not in neutral_loss_choices:
      continue
    peak_type = ions[i] * n_charges * n_neutral_losses + \
        neutral_loss_choices.index(neutral_losses[i]) * n_charges + \
        ion_charges[i] - 1
    y[positions[i]-1, peak_type] += intensities[i]
  y = y/np.sum(y)
  y = np.concatenate([y[:, :n_neutral_losses*n_charges],
                      np.flip(y[:, n_neutral_losses*n_charges:], 0)], 1)
  ys.append(y)
print(ys[0][0])
all_input = [(_X, _X_meta, _y) for _X, _X_meta, _y in zip(Xs, X_metas, ys)]




Xs = []
ys = []
X_metas = []

0  ['glycopeptide',
1 'charge',
2 'precursor_mass',
3 'retention_time',
4 'mzs',
5 'intensities',
6 'ions',
7 'positions',
8 'ion_charges']


# Read the path for input folder (CSV), then process the csv files into X, X_meta, and y,
# which are used in the deep learning model.
def read_csv(input_path):
    
    Read the CSV files from the folder of "input_path", then convert them into the formats for the deep learning model.
    :param input_path: A string for the input folder, such as "/data/Training-2021-D-Va--Human/CSV".
    :return:
    
    csv_path = Path(input_path)


    # sequence length
    seq_len = max_peptide_length + max_glycan_length
    for line in np.array(df):
        # Data for X: glycopeptide one-hot encoding (50, 26)
        Xs.append(np.array(eval(line[0])))
        # Data for X_meta: charge (1)
        X_metas.append(line[1])
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
        # The max position might equal to the "seq_len", so no need to "-1"
        # Use "b, y, Y" to represent ions, and times the maximum charges, then the space is 3*6 = 18
        y = np.zeros((seq_len, max_num_ions * max_num_charges))
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
                ion_type = (ion_charges[i] - 1) * max_num_ions + 1
            # Add the intensity to the corresponding place.
            y[position, ion_type] += intensities[i]
        # Normalization by the summary of the total intensities.
        y = y/np.sum(y)
        # Concatenate into y(50, 18)
        #y = np.concatenate([y[:, :n_neutral_losses*n_charges],
        #                  np.flip(y[:, n_neutral_losses*n_charges:], 0)], 1)
        ys.append(y)
    #print(len(ys[0]))
    #for i in range(len(ys[0])):
    #    print(ys[0][i])

    all_input = [(_X, _X_meta, _y) for _X, _X_meta, _y in zip(Xs, X_metas, ys)]
"""
# Change input_dim=24 to 26
net = TestModel(input_dim=26,
                n_tasks=max_num_ions * max_num_charges,
                embedding_dim=256,
                hidden_dim_lstm=128,
                hidden_dim_attention=32,
                n_lstm_layers=2,
                n_attention_heads=8,
                gpu=opt.gpu)
trainer = Trainer(net, opt)

#model_saved = "/home/zliang/orbitrap_data/data/Glycan-DeNovo-Test/MSP/model_bkp.pth"
#model_saved = "/home/zliang_venn_bio/spectranet_data/model_bkp.pth"
model_saved = "/data/saved_models/model_07_Mar.pth"
#trainer.load(model_saved)

#pickle_train_file ="/home/zliang_venn_bio/spectranet_data/train-66999_flipped.pkl"
#pickle_test_file ="/home/zliang_venn_bio/spectranet_data/test-9920_flipped.pkl"


pickle_train_path = "/data/Testing-01-Different-HCD-energies/Test01/PKL/"
pickle_train_file ="/data/Testing-01-Different-HCD-energies/Test01/PKL/210413-Energy-Test-B-15-20-35-40-1_HCDFT.pkl"
pickle_test_file ="/data/Testing-01-Different-HCD-energies/Test01/PKL/210413-Energy-Test-B-15-20-35-40-2_HCDFT.pkl"


train_input = []

pickle_train_folder = Path(pickle_train_path)

# Iterate all the files in the pickle folder
for pickle_train_file in pickle_train_folder.iterdir():
    train_pickle = pickle.load(open(pickle_train_file, 'rb'))
    # add each pickle file into the list, should use "extend", not "append".
    train_input.extend(train_pickle)


#train_input = pickle.load(open(pickle_train_file, 'rb'))
test_input = pickle.load(open(pickle_test_file, 'rb'))


"""
pd.set_option('display.width', None)
pd.set_option('display.max_rows',None)
pd.set_option('display.max_colwidth',None)

inf = str(test_input)
ft = open('test-9920_flipped_input.csv', 'w')
ft.write(inf)


df_train = pd.DataFrame(train_input)
df_train.to_csv(r'train-66999_flipped.csv')

df_test = pd.DataFrame(test_input)
df_test.to_csv(r'test-9920_flipped.csv')
"""

print("Training samples %d" % len(train_input))
print("Test samples %d" % len(test_input))

for i in range(8):
    start_fold_time = time.time()
    print("start fold %d" % i)
    if i > 0:
        evaluate(train_input, trainer.predict(train_input))
        evaluate(test_input, trainer.predict(test_input))
    end_fold_time = time.time()
    fold_time = end_fold_time - start_fold_time

    print("Time for evaluation:", fold_time)

    start_train_time = time.time()
    trainer.train(train_input, n_epochs=5)
    end_train_time = time.time()

    train_time = end_train_time - start_train_time
    print("Time for training:", train_time)

    if not os.path.exists('/data/saved_models'):
        os.mkdir('/data/saved_models')

    trainer.save(os.path.join('/data/saved_models', 'model-%d.pth' % i))

    end_save_time = time.time()
    save_time = end_save_time - end_train_time

    print("Time for saving:", save_time)
