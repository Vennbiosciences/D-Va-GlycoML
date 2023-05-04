#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 345 + 76 training files are from Training-all21-N-Human345_Cyno76-Full_Semi and Training-02-Cyno-76 .
And use transfer learning for Training-all21-O-Human345

The testing files are from Testing-07-O-2205-CRC .

The batch size is 8, learning rates are 1.

Created on 26 May 2022.
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
test_1_mean = [0.924112, 0.931932, 0.927047, 0.934762, 0.929187, 0.929447, 0.933965, 0.933270, 0.933781, 0.931827,
               0.934088, 0.935533, 0.934678, 0.932789, 0.933810, 0.933529, 0.933081, 0.931772, 0.934077, 0.929860,
               0.931096, 0.934260, 0.933077, 0.931546, 0.933792, 0.933137, 0.933661, 0.934466, 0.931927, 0.934041,
               0.934835, 0.933353, 0.933435, 0.932139, 0.932422, 0.931758, 0.934826, 0.933809, 0.933174, 0.932308,
               0.933062, 0.932660, 0.931327, 0.934355, 0.933102, 0.932867, 0.932518, 0.932693, 0.931883, 0.933212,
               ]
test_1_median = [0.965190, 0.971675, 0.966360, 0.974305, 0.965984, 0.964894, 0.970797, 0.972241, 0.972778, 0.969257,
                 0.974639, 0.973298, 0.971968, 0.973359, 0.973603, 0.972323, 0.971832, 0.971014, 0.972143, 0.968041,
                 0.971351, 0.970544, 0.970722, 0.972014, 0.973689, 0.971109, 0.971039, 0.973772, 0.969733, 0.972869,
                 0.970861, 0.972701, 0.971632, 0.968945, 0.971080, 0.967273, 0.973612, 0.973955, 0.971522, 0.972583,
                 0.970942, 0.972471, 0.971207, 0.972320, 0.972219, 0.972294, 0.971253, 0.971559, 0.967216, 0.971725,
                 ]

# The number of epochs for the training.
x = np.arange(1, 51)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(16, 8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1_mean, "r--.", label="Mean for Testing-07-O-2205-CRC on Transfer Learning Training-all21-O-Human345")
ax.plot(x, test_1_median, "g:+", label="Median for Testing-07-O-2205-CRC on Transfer Learning Training-all21-O-Human345")

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