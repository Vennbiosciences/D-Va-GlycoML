#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 345 + 76 training files are from Training-01-Human-345 and Training-02-Cyno-76 .
For Full and Semi.

The testing files are from Testing-03-2202-P4M .

The batch sizes are 32, 64, learning rates are 2, 4.

Created on 31 May 2022.
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
test_1_32_median = [0.958982, 0.963259, 0.965807, 0.967443, 0.966243, 0.965236, 0.967302, 0.967463, 0.968633, 0.968356,
                    0.968530, 0.968374, 0.967558, 0.968596, 0.968063, 0.968347, 0.967944, 0.967821, 0.969832, 0.968928,
                    0.969140, 0.968406, 0.968745, 0.968330, 0.968866, 0.968559, 0.968995, 0.968428, 0.968632, 0.969088,
                    0.969112, 0.969120, 0.968298, 0.968812, 0.968605, 0.968732, 0.968171, 0.968552, 0.969351, 0.969000,
                    0.969367, 0.969453, 0.969698, 0.968810, 0.968240, 0.968787, 0.969522, 0.969290, 0.968566, 0.968861,
                    ]
test_1_64_median = [0.960070, 0.964057, 0.963783, 0.966060, 0.967064, 0.966375, 0.967857, 0.968145, 0.968304, 0.968864,
                    0.967695, 0.968732, 0.968360, 0.968141, 0.966940, 0.967844, 0.968273, 0.967710, 0.968222, 0.968161,
                    0.968450, 0.969008, 0.968307, 0.968775, 0.968736, 0.969068, 0.969111, 0.968741, 0.969132, 0.968824,
                    0.968565, 0.968481, 0.967856, 0.967034, 0.967988, 0.968545, 0.968882, 0.967922, 0.968708, 0.968860,
                    0.968349, 0.968770, 0.968867, 0.969285, 0.967915, 0.968905, 0.968852, 0.968961, 0.969056, 0.969176,
                    ]
test_2_32_median = [0.965241, 0.966915, 0.968446, 0.966637, 0.968816, 0.968604, 0.968041, 0.968794, 0.969612, 0.969847,
                    0.969450, 0.968274, 0.969839, 0.969958, 0.968893, 0.969782, 0.968330, 0.969623, 0.968887, 0.968740,
                    0.969932, 0.970294, 0.970073, 0.969618, 0.969631, 0.970421, 0.970013, 0.969938, 0.970378, 0.970283,
                    0.969234, 0.969541, 0.969696, 0.969911, 0.969961, 0.970100, 0.970050, 0.970123, 0.969651, 0.970480,
                    0.969781, 0.969236, 0.970057, 0.970179, 0.969769, 0.970136, 0.969709, 0.970218, 0.970094, 0.969968,
                    ]
test_2_64_median = [0.962889, 0.965888, 0.967731, 0.968651, 0.967687, 0.969037, 0.968989, 0.968648, 0.968705, 0.969347,
                    0.969402, 0.967867, 0.968562, 0.970097, 0.969297, 0.969718, 0.969716, 0.970442, 0.969984, 0.969305,
                    0.969985, 0.969115, 0.970216, 0.969597, 0.968902, 0.970247, 0.968812, 0.969794, 0.969793, 0.970175,
                    0.968988, 0.969631, 0.969857, 0.969755, 0.970252, 0.968703, 0.970169, 0.970377, 0.970082, 0.969069,
                    0.969886, 0.970147, 0.970641, 0.969255, 0.970229, 0.969369, 0.969941, 0.970278, 0.969707, 0.969026,
                    ]
test_3_32_median = [0.963093, 0.966003, 0.967805, 0.967821, 0.967865, 0.968319, 0.969378, 0.969965, 0.969173, 0.969634,
                    0.968481, 0.969993, 0.968807, 0.969449, 0.969438, 0.970029, 0.969834, 0.969761, 0.968834, 0.969542,
                    0.968894, 0.970561, 0.969015, 0.970343, 0.968826, 0.969621, 0.969929, 0.969820, 0.969408, 0.969825,
                    0.969685, 0.970282, 0.970051, 0.969729, 0.968837, 0.969637, 0.969745, 0.969750, 0.970163, 0.970113,
                    0.970444, 0.969891, 0.969364, 0.970230, 0.969909, 0.970042, 0.968604, 0.969520, 0.969591, 0.969895,
               ]
test_3_64_median = [0.962030, 0.967776, 0.967101, 0.966822, 0.968295, 0.968577, 0.969008, 0.968778, 0.969341, 0.969041,
                    0.968495, 0.969488, 0.969012, 0.969001, 0.969203, 0.969649, 0.968477, 0.968089, 0.970413, 0.968756,
                    0.969249, 0.969194, 0.969603, 0.969329, 0.969466, 0.969384, 0.969712, 0.969811, 0.968484, 0.969240,
                    0.969214, 0.970105, 0.969833, 0.970191, 0.970359, 0.969384, 0.969880, 0.970208, 0.969995, 0.969284,
                    0.969854, 0.969714, 0.968943, 0.970063, 0.969496, 0.969771, 0.970100, 0.969032, 0.969328, 0.968868,
                    ]

# The number of epochs for the training.
x = np.arange(1, 51)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(16, 8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1_32_median, "r--.", label="Median for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76), BS=32, LR=2")
ax.plot(x, test_1_64_median, "r-x", label="Median for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76), BS=64, LR=4")
ax.plot(x, test_2_32_median, "b-.^", label="Median for Testing-03-2202-P4M on on Training-2021-Full(Human-345+Cyno-76), BS=32, LR=2")
ax.plot(x, test_2_64_median, "b:+", label="Median for Testing-03-2202-P4M on Training-2021-Full(Human-345+Cyno-76), BS=64, LR=4")
ax.plot(x, test_3_32_median, "y--+", label="Median for Testing-03-2202-P4M on Training-2021-Full-Semi(Human-345+Cyno-76), BS=32, LR=2")
ax.plot(x, test_3_64_median, "y-.,", label="Median for Testing-03-2202-P4M on Training-2021-Full-Semi(Human-345+Cyno-76) BS=64, LR=4")

# Set the scales for x-axis and y-axis.
ax.set_xlim([1, 50.1])
ax.set_ylim([0.90, 1])

# Set the scales for the display of x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 50
ax.set_xticks(np.linspace(1, 50, 50))
ax.set_yticks(np.linspace(0.90, 1, 3))

# Set the labels for the scales
ax.set_xticklabels(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
                    "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34",
                    "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"],
                   fontproperties="monospace", fontsize=12, rotation=10)
ax.set_yticklabels(["0.90", "0.95", "1.0"])

# Set the decorations for the scales and labels
ax.tick_params(left=False,pad=8,direction="in",length=2,width=3,color="b",labelsize=12)
ax.tick_params("x",labelrotation=10)

# Set the labels for x and y
ax.set_xlabel("Epochs", fontsize = 16)
ax.set_ylabel("Cos Similarity", fontsize = 16)

# Set the title
ax.set_title("Figure for Epochs and Cosine Similarities, Batch Size(32/64), Learning Rate(2/4)", fontsize=18,backgroundcolor='#F8CECC',
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