#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
################################################################################
This script create models based on Pytorch.

Created on 8 November 2021.
################################################################################
"""
__author__ = 'ZLiang'

from torch import nn
import torch as t
import torch.nn.functional as F
import numpy as np
from biLSTM import BiLSTM, MultiheadAttention, CombinedConv1D


#def CELoss(labels, outs, batch_size=8):
#def CELoss(labels, outs, batch_size=16):
def CELoss(labels, outs, batch_size=32):
#def CELoss(labels, outs, batch_size=64):
#def CELoss(labels, outs, batch_size=128):
#def CELoss(labels, outs, batch_size=256):
#def CELoss(labels, outs, batch_size=512):
    logits = outs.transpose(0, 1).contiguous().view(batch_size, -1)
    labels = labels.transpose(0, 1).contiguous().view((batch_size, -1))
    log_prob = F.log_softmax(logits, 1)
    return -(labels * log_prob).sum()


class TestModel(nn.Module):
    def __init__(self,
                 input_dim=20,
                 n_tasks=6,
                 embedding_dim=256,
                 hidden_dim_lstm=128,
                 hidden_dim_attention=32,
                 n_lstm_layers=2,
                 n_attention_heads=8,
                 gpu=True):
        """ Default reference-free model with LSTM and attention layers

        input_dim: int, optional
            size of input feature dimension
        hidden_dim_lstm: int, optional
            size of hidden LSTM layer
        hidden_dim_attention: int, optional
            size of each head in the attention layer
        n_lstm_layers: int, optional
            number of LSTM layers
        n_attention_heads: int, optional
            number of attention layer heads
        gpu: bool, optional
            if the model is run on GPU
        random_init: bool, optional
            if the initialize the LSTM hidden state randomly, for debug use
        """
        super(TestModel, self).__init__()
        self.input_dim = input_dim
        self.n_tasks = n_tasks
        self.embedding_dim = embedding_dim
        self.hidden_dim_lstm = hidden_dim_lstm
        self.hidden_dim_attention = hidden_dim_attention
        self.n_lstm_layers = n_lstm_layers
        self.n_attention_heads = n_attention_heads
        self.gpu = gpu

        self.embedding_module = t.nn.Embedding(self.input_dim, self.embedding_dim)
        # Add 1 for meta data,
        self.lstm_module = BiLSTM(input_dim=self.embedding_dim + 1,
                                  hidden_dim=self.hidden_dim_lstm,
                                  n_layers=self.n_lstm_layers,
                                  gpu=self.gpu)
        #self.lstm_module = BiLSTM(input_dim=self.embedding_dim + 2,
        #                          hidden_dim=self.hidden_dim_lstm,
        #                          n_layers=self.n_lstm_layers,
        #                          gpu=self.gpu)
        self.att_module = MultiheadAttention(Q_dim=self.hidden_dim_lstm,
                                             V_dim=self.hidden_dim_lstm,
                                             head_dim=self.hidden_dim_attention,
                                             n_heads=self.n_attention_heads)
        #    self.conv_module = CombinedConv1D(input_dim=self.hidden_dim_lstm,
        #                                      conv_dim=self.hidden_dim_attention)
        self.fc = nn.Sequential(
            nn.Linear(self.hidden_dim_lstm * 2, 64),
            nn.ReLU(True),
            nn.Linear(64, self.n_tasks))

    def forward(self, sequence, metas, batch_size=1):
        """
        sequence: torch Tensor
            input batch of shape seq_len * batch_size * input_dim
        metas: precursor charge information for the encoding
        batch_size: int, optional
            batch size
        """
        # sequence of shape: seq_len * batch_size
        embedded_inputs = t.matmul(sequence, self.embedding_module.weight)
        seq_len = embedded_inputs.shape[0]
        # *metas.shape equals to metas.shape[0]
        #metas = metas.view(1, *metas.shape).repeat(seq_len, 1, 1)
        metas = metas.view(1, *metas.shape).repeat(seq_len, 1, 1)

        lstm_inputs = t.cat([embedded_inputs, metas], 2)
        lstm_outs = self.lstm_module(lstm_inputs, batch_size=batch_size)
        attention_outs = self.att_module(sequence=lstm_outs)
        # conv_outs = self.conv_module(lstm_outs)
        intervals = t.cat([attention_outs[1:], attention_outs[:-1]], 2)
        # outs of shape: (seq_len - 1) * batch_size * n_tasks
        return self.fc(intervals)

    def predict(self, sequence, metas, batch_size=1, gpu=True):
        """
        sequence: torch Tensor
            input batch of shape seq_len * batch_size * input_dim
        batch_size: int, optional
            batch size
        gpu: bool, optional
            if run on GPU
        """
        assert batch_size == 1
        output = self.forward(sequence, metas, 1)[:, 0, :]
        output_shape = output.shape
        logits = output.view(-1)
        log_prob = F.log_softmax(logits, 0)
        prob = t.exp(log_prob).view(output_shape)
        if gpu:
            prob = prob.cpu()
        prob = prob.data.numpy()
        return prob
