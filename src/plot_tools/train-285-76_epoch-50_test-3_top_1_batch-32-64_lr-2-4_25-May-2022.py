#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 285 + 76 training files are from Training-01-Human-285 and Training-02-Cyno-76 .

The testing files are from Testing-03-2202-P4M .

The batch size is 32, learning rates are 2.
The batch size is 64, learning rates are 4.

Created on 25 May 2022.
###################################################################################################################
"""
__author__ = 'ZLiang'

import matplotlib.pyplot as plt
import numpy as np

# Cosine similarities for different test cases. Totally 50 scores for each test.
# These are used for the y-axis.
# The values are from the epochs of
#   1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
#   27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50.
test_1_mean = [0.931371, 0.935344, 0.938644, 0.940911, 0.938898, 0.938229, 0.940869, 0.940547, 0.941919, 0.942604,
               0.941553, 0.941499, 0.940713, 0.942472, 0.940594, 0.942066, 0.940894, 0.940625, 0.943034, 0.942443,
               0.942420, 0.941801, 0.942640, 0.941968, 0.942130, 0.941952, 0.942593, 0.942453, 0.941989, 0.942294,
               0.942566, 0.942991, 0.941863, 0.942318, 0.942521, 0.942276, 0.941522, 0.942421, 0.942528, 0.942779,
               0.942746, 0.942916, 0.943084, 0.941984, 0.941974, 0.942720, 0.942437, 0.942749, 0.941874, 0.943232,
               ]
test_1_median = [0.958982, 0.963259, 0.965807, 0.967443, 0.966243, 0.965236, 0.967302, 0.967463, 0.968633, 0.968356,
                 0.968530, 0.968374, 0.967558, 0.968596, 0.968063, 0.968347, 0.967944, 0.967821, 0.969832, 0.968928,
                 0.969140, 0.968406, 0.968745, 0.968330, 0.968866, 0.968559, 0.968995, 0.968428, 0.968632, 0.969088,
                 0.969112, 0.969120, 0.968298, 0.968812, 0.968605, 0.968732, 0.968171, 0.968552, 0.969351, 0.969000,
                 0.969367, 0.969453, 0.969698, 0.968810, 0.968240, 0.968787, 0.969522, 0.969290, 0.968566, 0.968861,
                 ]
test_2_mean = [0.932984, 0.937450, 0.936825, 0.939107, 0.940033, 0.938686, 0.940594, 0.940474, 0.942007, 0.942482,
               0.941806, 0.941658, 0.941204, 0.942107, 0.940594, 0.940820, 0.941329, 0.940213, 0.941564, 0.941317,
               0.941233, 0.942423, 0.941341, 0.942392, 0.942575, 0.942361, 0.942091, 0.942192, 0.942374, 0.942455,
               0.941989, 0.942426, 0.941756, 0.941088, 0.941296, 0.942227, 0.942440, 0.941517, 0.942561, 0.942624,
               0.941731, 0.941710, 0.942436, 0.942819, 0.942112, 0.942364, 0.942640, 0.942014, 0.942084, 0.942621,
               ]
test_2_median = [0.960070, 0.964057, 0.963783, 0.966060, 0.967064, 0.966375, 0.967857, 0.968145, 0.968304, 0.968864,
                 0.967695, 0.968732, 0.968360, 0.968141, 0.966940, 0.967844, 0.968273, 0.967710, 0.968222, 0.968161,
                 0.968450, 0.969008, 0.968307, 0.968775, 0.968736, 0.969068, 0.969111, 0.968741, 0.969132, 0.968824,
                 0.968565, 0.968481, 0.967856, 0.967034, 0.967988, 0.968545, 0.968882, 0.967922, 0.968708, 0.968860,
                 0.968349, 0.968770, 0.968867, 0.969285, 0.967915, 0.968905, 0.968852, 0.968961, 0.969056, 0.969176,
                 ]

# The number of epochs for the training.
x = np.arange(1, 51)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(16, 8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1_mean, "r--.", label="Mean for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76), Batch=32, LR=0.0002")
ax.plot(x, test_1_median, "r-x", label="Median for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76), Batch=32, LR=0.0002")
ax.plot(x, test_2_mean, "b-.^", label="Mean for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76, Batch=64, LR=0.0004")
ax.plot(x, test_2_median, "b:+", label="Median for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76), Batch=64, LR=0.0004")

# Set the scales for x-axis and y-axis.
ax.set_xlim([1, 50.1])
ax.set_ylim([0.85, 1])

# Set the scales for the display of x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 50
ax.set_xticks(np.linspace(1, 50, 50))
ax.set_yticks(np.linspace(0.85, 1, 4))

# Set the labels for the scales
ax.set_xticklabels(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
                    "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34",
                    "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"],
                   fontproperties="monospace", fontsize=12, rotation=10)
ax.set_yticklabels(["0.85", "0.90", "0.95", "1.0"])

# Set the decorations for the scales and labels
ax.tick_params(left=False,pad=8,direction="in",length=2,width=3,color="b",labelsize=12)
ax.tick_params("x",labelrotation=10)

# Set the labels for x and y
ax.set_xlabel("Epochs", fontsize = 16)
ax.set_ylabel("Cos Similarity", fontsize = 16)

# Set the title
ax.set_title("Figure for Epochs and Cosine Similarities", fontsize=18,backgroundcolor='#DAE8FC',
             fontweight='bold',color='black',verticalalignment="baseline")

# Set the spines
ax.spines["left"].set_color("darkblue")
ax.spines["bottom"].set_linewidth(2)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

# Set the annotations
#ax.annotate(s="Max:0.960398",xy=(17,0.94),xytext=(17.2,0.87),arrowprops=dict(facecolor="g",shrink=0.05,headwidth=12,
#                                                                          headlength=6, width=4),fontsize=12)
#ax.annotate(s="Min:0.421085",xy=(1,0.4),xytext=(1.2,0.33),arrowprops=dict(facecolor="b",shrink=0.05,headwidth=12,
#                                                                          headlength=6, width=4),fontsize=12)

# Set the legend
ax.legend(loc=3,labelspacing=1,handlelength=3,fontsize=14,shadow=True)

plt.show()