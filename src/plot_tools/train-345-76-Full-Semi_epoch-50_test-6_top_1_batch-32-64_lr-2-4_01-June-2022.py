#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 345 + 76 training files are from Training-01-Human-345 and Training-02-Cyno-76 .
For Full and Semi.

The testing files are from Testing-06-N-2205-CRC-Semi .

The batch sizes are 32, 64, learning rates are 2, 4.

Created on 01 June 2022.
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
test_1_32_median = [0.907395, 0.904357, 0.909728, 0.913999, 0.904275, 0.905619, 0.907779, 0.904356, 0.897990, 0.906158,
                    0.907275, 0.909065, 0.900799, 0.907013, 0.904072, 0.903542, 0.903005, 0.901237, 0.907239, 0.907990,
                    0.903472, 0.906898, 0.903907, 0.902699, 0.907757, 0.905430, 0.902774, 0.906923, 0.901842, 0.906021,
                    0.908747, 0.906159, 0.906998, 0.903984, 0.907908, 0.904194, 0.902380, 0.906812, 0.904021, 0.909882,
                    0.905634, 0.906517, 0.906671, 0.905791, 0.904424, 0.903275, 0.908346, 0.907146, 0.903932, 0.906313,
                    ]
test_1_64_median = [0.897728, 0.905283, 0.901699, 0.903750, 0.901677, 0.902819, 0.901035, 0.899776, 0.897635, 0.899989,
                    0.899384, 0.897430, 0.901708, 0.900732, 0.896018, 0.896709, 0.894747, 0.894247, 0.900336, 0.898632,
                    0.896212, 0.898330, 0.897883, 0.895781, 0.894830, 0.892981, 0.895267, 0.895786, 0.895433, 0.895962,
                    0.896675, 0.894669, 0.897523, 0.898806, 0.892839, 0.899466, 0.895351, 0.890706, 0.896838, 0.896814,
                    0.895284, 0.889758, 0.892580, 0.895668, 0.892161, 0.895579, 0.895546, 0.894565, 0.896677, 0.892183,
                    ]
test_2_32_median = [0.889872, 0.892153, 0.900473, 0.902292, 0.907023, 0.901345, 0.902095, 0.899022, 0.904394, 0.905778,
                    0.903621, 0.898167, 0.900843, 0.902430, 0.900433, 0.903269, 0.898360, 0.903411, 0.903861, 0.902096,
                    0.906278, 0.905116, 0.905791, 0.896868, 0.902448, 0.906292, 0.904861, 0.900395, 0.905151, 0.902915,
                    0.900008, 0.903797, 0.903757, 0.902887, 0.903429, 0.903534, 0.898832, 0.903040, 0.901095, 0.904402,
                    0.901966, 0.898776, 0.891617, 0.905075, 0.900695, 0.901336, 0.902037, 0.899205, 0.899179, 0.898978,
                    ]
test_2_64_median = [0.902516, 0.903948, 0.910142, 0.909787, 0.905175, 0.905890, 0.906979, 0.901116, 0.901454, 0.904700,
                    0.900509, 0.902681, 0.905110, 0.907744, 0.907915, 0.905404, 0.906496, 0.905816, 0.906510, 0.907084,
                    0.906682, 0.906605, 0.907800, 0.907532, 0.908362, 0.907946, 0.905151, 0.906685, 0.905597, 0.905763,
                    0.905603, 0.904443, 0.904011, 0.905940, 0.908444, 0.903161, 0.909731, 0.909712, 0.908985, 0.907342,
                    0.908673, 0.908451, 0.906355, 0.907383, 0.908271, 0.905001, 0.907218, 0.904831, 0.906290, 0.902910,
                    ]
test_3_32_median = [0.949816, 0.957600, 0.958447, 0.960120, 0.958150, 0.959109, 0.958647, 0.961661, 0.957812, 0.961265,
                    0.958944, 0.961579, 0.961769, 0.962882, 0.960619, 0.962079, 0.961141, 0.961269, 0.956423, 0.962752,
                    0.961340, 0.960493, 0.961079, 0.961401, 0.960241, 0.961618, 0.961614, 0.960467, 0.960325, 0.960831,
                    0.961808, 0.962153, 0.962789, 0.960214, 0.961482, 0.961459, 0.960482, 0.961896, 0.959757, 0.961945,
                    0.961588, 0.961112, 0.960994, 0.962114, 0.961485, 0.962451, 0.959263, 0.960653, 0.960555, 0.961506,
                    ]
test_3_64_median = [0.952039, 0.957417, 0.960140, 0.956097, 0.961017, 0.959590, 0.958992, 0.958642, 0.958923, 0.960858,
                    0.959715, 0.960794, 0.957794, 0.959544, 0.960939, 0.960965, 0.957740, 0.957448, 0.961859, 0.958926,
                    0.961365, 0.958337, 0.958819, 0.960673, 0.961141, 0.960007, 0.960925, 0.960112, 0.958550, 0.960317,
                    0.962945, 0.962384, 0.961008, 0.959710, 0.960944, 0.960283, 0.960262, 0.962108, 0.960905, 0.960384,
                    0.962854, 0.961123, 0.960991, 0.963496, 0.959815, 0.961642, 0.960239, 0.961085, 0.960843, 0.959989,
                    ]

# The number of epochs for the training.
x = np.arange(1, 51)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(16, 8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1_32_median, "r--.", label="Median for Testing-06-N-2205-CRC-Semi on Training-02(Human-285+Cyno-76), BS=32, LR=2")
ax.plot(x, test_1_64_median, "r-x", label="Median for Testing-06-N-2205-CRC-Semi on Training-02(Human-285+Cyno-76), BS=64, LR=4")
ax.plot(x, test_2_32_median, "b-.^", label="Median for Testing-06-N-2205-CRC-Semi on on Training-2021-Full(Human-345+Cyno-76), BS=32, LR=2")
ax.plot(x, test_2_64_median, "b:+", label="Median for Testing-06-N-2205-CRC-Semi on Training-2021-Full(Human-345+Cyno-76), BS=64, LR=4")
ax.plot(x, test_3_32_median, "y--+", label="Median for Testing-06-N-2205-CRC-Semi on Training-2021-Full-Semi(Human-345+Cyno-76), BS=32, LR=2")
ax.plot(x, test_3_64_median, "y-.,", label="Median for Testing-06-N-2205-CRC-Semi on Training-2021-Full-Semi(Human-345+Cyno-76), BS=64, LR=4")

# Set the scales for x-axis and y-axis.
ax.set_xlim([1, 50.1])
ax.set_ylim([0.80, 1])

# Set the scales for the display of x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 50
ax.set_xticks(np.linspace(1, 50, 50))
ax.set_yticks(np.linspace(0.80, 1, 5))

# Set the labels for the scales
ax.set_xticklabels(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
                    "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34",
                    "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"],
                   fontproperties="monospace", fontsize=12, rotation=10)
ax.set_yticklabels(["0.80", "0.85", "0.90", "0.95", "1.0"])

# Set the decorations for the scales and labels
ax.tick_params(left=False,pad=8,direction="in",length=2,width=3,color="b",labelsize=12)
ax.tick_params("x",labelrotation=10)

# Set the labels for x and y
ax.set_xlabel("Epochs", fontsize = 16)
ax.set_ylabel("Cos Similarity", fontsize = 16)

# Set the title
ax.set_title("Figure for Epochs and Cosine Similarities, Batch Size(32/64), Learning Rate(2/4)", fontsize=18,backgroundcolor='#E1D5E7',
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