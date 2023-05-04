#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 345 + 76 training files are from Training-01-Human-345 and Training-02-Cyno-76 .
For Full and Semi.

The testing files are from Testing-05-N-2205-CRC .

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
test_1_32_median = [0.947845, 0.952665, 0.955552, 0.956979, 0.954820, 0.953996, 0.956639, 0.955426, 0.957225, 0.955983,
                    0.957513, 0.958119, 0.956085, 0.957136, 0.957277, 0.957403, 0.957090, 0.956847, 0.958656, 0.957944,
                    0.957691, 0.957871, 0.957527, 0.957494, 0.958387, 0.957512, 0.957001, 0.957765, 0.957036, 0.958399,
                    0.958316, 0.957907, 0.957389, 0.957770, 0.957589, 0.957728, 0.956783, 0.957574, 0.957915, 0.957893,
                    0.957905, 0.958242, 0.958549, 0.958085, 0.957078, 0.957782, 0.957792, 0.958181, 0.957732, 0.958018,
                    ]
test_1_64_median = [0.948904, 0.952884, 0.952349, 0.954641, 0.955204, 0.955141, 0.956971, 0.956121, 0.956737, 0.957736,
                    0.956624, 0.957066, 0.958114, 0.957205, 0.956366, 0.957182, 0.956862, 0.957043, 0.957689, 0.957673,
                    0.957516, 0.958885, 0.957320, 0.957228, 0.957547, 0.957838, 0.957904, 0.957269, 0.958268, 0.957544,
                    0.956924, 0.957358, 0.957510, 0.956797, 0.956226, 0.957901, 0.957479, 0.956795, 0.957907, 0.958123,
                    0.957576, 0.957677, 0.957553, 0.957787, 0.956882, 0.957837, 0.958148, 0.958054, 0.958159, 0.957999,
                    ]
test_2_32_median = [0.952954, 0.953623, 0.955957, 0.955209, 0.956011, 0.955847, 0.955909, 0.956638, 0.957728, 0.957566,
                    0.957595, 0.955674, 0.957544, 0.957917, 0.957189, 0.957941, 0.956780, 0.957839, 0.957622, 0.956861,
                    0.958543, 0.958338, 0.958983, 0.957211, 0.957829, 0.958643, 0.958147, 0.957252, 0.958978, 0.957551,
                    0.956695, 0.957932, 0.958230, 0.958273, 0.957903, 0.958435, 0.957129, 0.958458, 0.957947, 0.958682,
                    0.958055, 0.957738, 0.957406, 0.958626, 0.957816, 0.958122, 0.958310, 0.957705, 0.958040, 0.958282,
                    ]
test_2_64_median = [0.949081, 0.952765, 0.955474, 0.955899, 0.955781, 0.957382, 0.956748, 0.956455, 0.956706, 0.957139,
                    0.957744, 0.956144, 0.956948, 0.957284, 0.957573, 0.957690, 0.958063, 0.958262, 0.957959, 0.957129,
                    0.957975, 0.957289, 0.958456, 0.957387, 0.956677, 0.958359, 0.957404, 0.957974, 0.957603, 0.958318,
                    0.957335, 0.957930, 0.957934, 0.957935, 0.958293, 0.956700, 0.958345, 0.959207, 0.958376, 0.957005,
                    0.958056, 0.958530, 0.957911, 0.957828, 0.958492, 0.957685, 0.958355, 0.958172, 0.957935, 0.957580,
                    ]
test_3_32_median = [0.950620, 0.953782, 0.955024, 0.955994, 0.955027, 0.955891, 0.957244, 0.957531, 0.956551, 0.956847,
                    0.956304, 0.958078, 0.957054, 0.957855, 0.957144, 0.957309, 0.957301, 0.957046, 0.956049, 0.957880,
                    0.957493, 0.958024, 0.957250, 0.958944, 0.957388, 0.958101, 0.957926, 0.957693, 0.957658, 0.957598,
                    0.957674, 0.957924, 0.958207, 0.957213, 0.957256, 0.958117, 0.957765, 0.957914, 0.958281, 0.958386,
                    0.957985, 0.957619, 0.957979, 0.958933, 0.958585, 0.958263, 0.956725, 0.958249, 0.957593, 0.958066,
                    ]
test_3_64_median = [0.948308, 0.955290, 0.954514, 0.954122, 0.956259, 0.956485, 0.957261, 0.956653, 0.956855, 0.957456,
                    0.956497, 0.956859, 0.956929, 0.957336, 0.956783, 0.957357, 0.955985, 0.956078, 0.958772, 0.956351,
                    0.957516, 0.957046, 0.957473, 0.957021, 0.957779, 0.958034, 0.957648, 0.957027, 0.956826, 0.956589,
                    0.957368, 0.957934, 0.957622, 0.957731, 0.958384, 0.957346, 0.957899, 0.958627, 0.957753, 0.957640,
                    0.958226, 0.958437, 0.957332, 0.958140, 0.957832, 0.957312, 0.957846, 0.957585, 0.957900, 0.957421,
                    ]

# The number of epochs for the training.
x = np.arange(1, 51)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(16, 8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1_32_median, "r--.", label="Median for Testing-05-N-2205-CRC on Training-02(Human-285+Cyno-76), BS=32, LR=2")
ax.plot(x, test_1_64_median, "r-x", label="Median for Testing-05-N-2205-CRC on Training-02(Human-285+Cyno-76), BS=64, LR=4")
ax.plot(x, test_2_32_median, "b-.^", label="Median for Testing-05-N-2205-CRC on on Training-2021-Full(Human-345+Cyno-76), BS=32, LR=2")
ax.plot(x, test_2_64_median, "b:+", label="Median for Testing-05-N-2205-CRC on Training-2021-Full(Human-345+Cyno-76), BS=64, LR=4")
ax.plot(x, test_3_32_median, "y--+", label="Median for Testing-05-N-2205-CRC on Training-2021-Full-Semi(Human-345+Cyno-76), BS=32, LR=2")
ax.plot(x, test_3_64_median, "y-.,", label="Median for Testing-05-N-2205-CRC on Training-2021-Full-Semi(Human-345+Cyno-76), BS=64, LR=4")

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
ax.set_title("Figure for Epochs and Cosine Similarities, Batch Size(32/64), Learning Rate(2/4)", fontsize=18,backgroundcolor='#FFE6CC',
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