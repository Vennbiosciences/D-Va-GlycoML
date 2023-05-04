#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 345 + 76 training files are from Training-01-Human-345 and Training-02-Cyno-76 .
For Full and Semi.

The testing files are from Testing-03-2202-P4M .

The batch size is 64, learning rates is 4.

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
test_1_mean = [0.932984, 0.937450, 0.936825, 0.939107, 0.940033, 0.938686, 0.940594, 0.940474, 0.942007, 0.942482,
               0.941806, 0.941658, 0.941204, 0.942107, 0.940594, 0.940820, 0.941329, 0.940213, 0.941564, 0.941317,
               0.941233, 0.942423, 0.941341, 0.942392, 0.942575, 0.942361, 0.942091, 0.942192, 0.942374, 0.942455,
               0.941989, 0.942426, 0.941756, 0.941088, 0.941296, 0.942227, 0.942440, 0.941517, 0.942561, 0.942624,
               0.941731, 0.941710, 0.942436, 0.942819, 0.942112, 0.942364, 0.942640, 0.942014, 0.942084, 0.942621,
               ]
test_1_median = [0.960070, 0.964057, 0.963783, 0.966060, 0.967064, 0.966375, 0.967857, 0.968145, 0.968304, 0.968864,
                 0.967695, 0.968732, 0.968360, 0.968141, 0.966940, 0.967844, 0.968273, 0.967710, 0.968222, 0.968161,
                 0.968450, 0.969008, 0.968307, 0.968775, 0.968736, 0.969068, 0.969111, 0.968741, 0.969132, 0.968824,
                 0.968565, 0.968481, 0.967856, 0.967034, 0.967988, 0.968545, 0.968882, 0.967922, 0.968708, 0.968860,
                 0.968349, 0.968770, 0.968867, 0.969285, 0.967915, 0.968905, 0.968852, 0.968961, 0.969056, 0.969176,
                 ]
test_2_mean = [0.936067, 0.938631, 0.941064, 0.942160, 0.940842, 0.942009, 0.942938, 0.942072, 0.942470, 0.942626,
               0.943057, 0.941419, 0.942406, 0.944140, 0.942829, 0.943575, 0.942909, 0.944181, 0.943686, 0.943859,
               0.944136, 0.943242, 0.944594, 0.944313, 0.943458, 0.944272, 0.943315, 0.943990, 0.944318, 0.944282,
               0.943574, 0.943773, 0.944216, 0.943105, 0.944763, 0.942367, 0.944435, 0.944562, 0.944654, 0.943014,
               0.943607, 0.944114, 0.944975, 0.943517, 0.944177, 0.944175, 0.944131, 0.944034, 0.944501, 0.942644,
               ]
test_2_median = [0.962889, 0.965888, 0.967731, 0.968651, 0.967687, 0.969037, 0.968989, 0.968648, 0.968705, 0.969347,
                 0.969402, 0.967867, 0.968562, 0.970097, 0.969297, 0.969718, 0.969716, 0.970442, 0.969984, 0.969305,
                 0.969985, 0.969115, 0.970216, 0.969597, 0.968902, 0.970247, 0.968812, 0.969794, 0.969793, 0.970175,
                 0.968988, 0.969631, 0.969857, 0.969755, 0.970252, 0.968703, 0.970169, 0.970377, 0.970082, 0.969069,
                 0.969886, 0.970147, 0.970641, 0.969255, 0.970229, 0.969369, 0.969941, 0.970278, 0.969707, 0.969026,
                 ]
test_3_mean = [0.935344, 0.940521, 0.940336, 0.939621, 0.941976, 0.942741, 0.941880, 0.941760, 0.942740, 0.942813,
               0.942249, 0.943203, 0.942829, 0.942938, 0.943691, 0.942867, 0.941504, 0.941678, 0.943225, 0.943187,
               0.942863, 0.943573, 0.943458, 0.943321, 0.943898, 0.943362, 0.943536, 0.943747, 0.941967, 0.943148,
               0.943154, 0.944081, 0.944122, 0.943133, 0.944109, 0.943544, 0.944028, 0.944134, 0.944269, 0.943438,
               0.944193, 0.943731, 0.943662, 0.944738, 0.943230, 0.944320, 0.944050, 0.943078, 0.942789, 0.942971,
               ]
test_3_median = [0.962030, 0.967776, 0.967101, 0.966822, 0.968295, 0.968577, 0.969008, 0.968778, 0.969341, 0.969041,
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
ax.plot(x, test_1_mean, "r--.", label="Mean for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76)")
ax.plot(x, test_1_median, "r-x", label="Median for Testing-03-2202-P4M on Training-02(Human-285+Cyno-76)")
ax.plot(x, test_2_mean, "b-.^", label="Mean for Testing-03-2202-P4M on on Training-2021-Full(Human-345+Cyno-76)")
ax.plot(x, test_2_median, "b:+", label="Median for Testing-03-2202-P4M on Training-2021-Full(Human-345+Cyno-76)")
ax.plot(x, test_3_mean, "y--+", label="Mean for Testing-03-2202-P4M on Training-2021-Full-Semi(Human-345+Cyno-76)")
ax.plot(x, test_3_median, "y-.,", label="Median for Testing-03-2202-P4M on Training-2021-Full-Semi(Human-345+Cyno-76)")

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
ax.set_title("Figure for Epochs and Cosine Similarities, Batch Size(64), Learning Rate(4)", fontsize=18,backgroundcolor='#F8CECC',
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