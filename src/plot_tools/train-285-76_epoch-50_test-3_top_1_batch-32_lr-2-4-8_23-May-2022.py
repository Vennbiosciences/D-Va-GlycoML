#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 285 + 76 training files are from Training-01-Human-285 and Training-02-Cyno-76 .

The testing files are from Testing-03-2202-P4M .

The batch size is 32, learning rates are 2, 4, 8.

Created on 23 May 2022.
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
test_2_mean = [0.931473, 0.937870, 0.938423, 0.939440, 0.939583, 0.940636, 0.940890, 0.940839, 0.941257, 0.940671,
               0.940692, 0.941430, 0.941660, 0.941147, 0.942217, 0.940917, 0.941267, 0.941283, 0.941444, 0.941381,
               0.941602, 0.941646, 0.942301, 0.942007, 0.941376, 0.942308, 0.942023, 0.941271, 0.942519, 0.941145,
               0.940780, 0.942524, 0.942756, 0.940195, 0.942261, 0.942083, 0.941801, 0.942581, 0.941842, 0.942157,
               0.941743, 0.941841, 0.942406, 0.940792, 0.942217, 0.941386, 0.941786, 0.942160, 0.940959, 0.941595,
               ]
test_2_median = [0.961798, 0.965290, 0.966489, 0.966243, 0.965702, 0.966666, 0.968620, 0.968375, 0.968039, 0.967472,
                 0.968113, 0.967965, 0.968256, 0.968626, 0.968251, 0.968706, 0.967747, 0.968420, 0.968191, 0.967892,
                 0.968352, 0.968044, 0.969043, 0.968071, 0.967867, 0.968870, 0.968266, 0.968343, 0.969165, 0.967472,
                 0.968079, 0.968560, 0.969074, 0.966931, 0.969231, 0.969055, 0.968214, 0.969100, 0.968589, 0.968756,
                 0.968666, 0.968323, 0.969047, 0.967791, 0.968800, 0.968387, 0.968211, 0.968288, 0.967901, 0.967809,
                 ]
test_3_mean = [0.936050, 0.936552, 0.938620, 0.938269, 0.940299, 0.940129, 0.938349, 0.939704, 0.940729, 0.941119,
               0.940519, 0.940450, 0.941081, 0.941214, 0.939998, 0.940836, 0.940660, 0.940972, 0.941088, 0.941898,
               0.941082, 0.941618, 0.940378, 0.941266, 0.941560, 0.940178, 0.941799, 0.939456, 0.940629, 0.941672,
               0.940564, 0.941747, 0.941462, 0.940776, 0.939401, 0.941371, 0.939547, 0.942387, 0.938610, 0.940781,
               0.939031, 0.940512, 0.941007, 0.941758, 0.940004, 0.938145, 0.940315, 0.941990, 0.941726, 0.940222,
               ]
test_3_median = [0.964611, 0.964284, 0.966489, 0.965420, 0.967054, 0.967514, 0.965648, 0.967183, 0.967516, 0.968061,
                 0.967514, 0.966450, 0.967614, 0.967992, 0.967258, 0.968243, 0.967605, 0.968618, 0.967435, 0.967984,
                 0.967614, 0.967380, 0.967025, 0.967255, 0.967633, 0.966398, 0.968156, 0.965927, 0.967663, 0.968031,
                 0.967787, 0.966317, 0.968020, 0.967590, 0.966094, 0.968285, 0.967425, 0.968917, 0.966112, 0.967437,
                 0.965521, 0.967291, 0.967089, 0.967429, 0.968018, 0.965639, 0.967977, 0.967705, 0.968314, 0.966459,
                 ]

# The number of epochs for the training.
x = np.arange(1, 51)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(16, 8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1_mean, "r--.", label="Mean for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76), LR=0.0002")
ax.plot(x, test_1_median, "r-x", label="Median for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76), LR=0.0002")
ax.plot(x, test_2_mean, "b-.^", label="Mean for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76, LR=0.0004")
ax.plot(x, test_2_median, "b:+", label="Median for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76), LR=0.0004")
ax.plot(x, test_3_mean, "y--+", label="Mean for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76, LR=0.0008")
ax.plot(x, test_3_median, "y-.,", label="Median for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76), LR=0.0008")

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
ax.set_title("Figure for Epochs and Cosine Similarities", fontsize=18,backgroundcolor='#FFF2CC',
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