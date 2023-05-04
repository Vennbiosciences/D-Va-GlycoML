#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 15 training files are from Energy-01-HCD-15-20-34-37-40 .
The 2 testing file are from Energy-01-HCD-15-20-34-37-40 with 2 different batch size (8 and 32).

Created on 18 April 2022.
###################################################################################################################
"""
__author__ = 'ZLiang'

import matplotlib.pyplot as plt
import numpy as np

# Cosine similarities for different test cases. Totally 34 scores for each test.
# These are used for the y-axis.
# The values are from the folds of
#   1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49,
#   52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 100
test_1_mean = [0.924370, 0.934690, 0.941953, 0.945356, 0.946473, 0.947832, 0.949238, 0.949706, 0.950814, 0.951501,
               0.951591, 0.952122, 0.952589, 0.952236, 0.952632, 0.952517, 0.953037, 0.953628, 0.953245, 0.953362,
               0.953520, 0.954029, 0.953786, 0.953944, 0.953531, 0.953898, 0.954066, 0.954378, 0.954429, 0.954155,
               0.953974, 0.954168, 0.953711, 0.954282, ]
test_1_median = [0.963244, 0.970500, 0.976539, 0.978251, 0.979320, 0.980713, 0.979787, 0.980811, 0.982511, 0.982191,
                 0.982257, 0.982328, 0.982549, 0.982819, 0.983529, 0.982722, 0.984188, 0.983632, 0.983864, 0.983869,
                 0.983656, 0.983801, 0.983186, 0.984468, 0.982545, 0.983873, 0.984597, 0.984414, 0.982611, 0.984495,
                 0.984998, 0.984456, 0.984699, 0.985530, ]
test_2_mean = [0.906497, 0.930046, 0.934685, 0.940447, 0.941678, 0.943164, 0.944249, 0.945662, 0.946849, 0.947162,
               0.948006, 0.948585, 0.949558, 0.949766, 0.950064, 0.950493, 0.950267, 0.950190, 0.951273, 0.951662,
               0.951730, 0.951904, 0.952168, 0.951894, 0.951936, 0.952704, 0.953261, 0.953143, 0.952336, 0.952828,
               0.953062, 0.953390, 0.953276, 0.953533, ]
test_2_median = [0.947792, 0.966229, 0.968565, 0.973638, 0.974975, 0.976592, 0.976935, 0.976845, 0.978284, 0.977891,
                 0.978531, 0.981048, 0.980581, 0.980989, 0.981515, 0.981616, 0.981804, 0.981979, 0.981936, 0.982450,
                 0.982224, 0.981515, 0.982667, 0.982805, 0.983837, 0.982686, 0.983851, 0.983352, 0.983186, 0.983982,
                 0.984352, 0.984302, 0.984057, 0.983349, ]

# The number of folds for the training, each fold contains 5 epochs.
# The folds range from 1 to 35, corresponding the epochs range from 5 to 500.
x = np.arange(1, 35)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(16, 8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1_mean, "b--.", label="Testing the Mean for Energy-01-HCD-15-20-34-37-40 with batch size 8")
ax.plot(x, test_1_median, "r-x", label="Testing the Median for Energy-01-HCD-15-20-34-37-40 with batch size 8")
ax.plot(x, test_2_mean, "y-.^", label="Testing the Mean for Energy-02-HCD-15-20-35-40 with batch size 32")
ax.plot(x, test_2_median, "g:+", label="Testing the Median for Energy-02-HCD-15-20-35-40 with batch size 32")

# Set the scales for x-axis and y-axis.
ax.set_xlim([1, 34.1])
ax.set_ylim([0.8, 1])

# Set the scales for the display of x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 34
ax.set_xticks(np.linspace(1, 34, 34))
ax.set_yticks(np.linspace(0.8, 1, 3))

# Set the labels for the scales
ax.set_xticklabels(["5", "20", "35", "50", "65", "80", "95", "110", "125", "140", "155", "170", "185", "200", "215",
                    "230", "245", "260", "275", "290", "305", "320", "335", "350", "365", "380", "395", "410", "425",
                    "440", "455", "470", "485", "500", ], fontproperties="monospace", fontsize=12, rotation=10)
ax.set_yticklabels(["0.8", "0.9", "1.0"])

# Set the decorations for the scales and labels
ax.tick_params(left=False,pad=8,direction="in",length=2,width=3,color="b",labelsize=12)
ax.tick_params("x",labelrotation=10)

# Set the labels for x and y
ax.set_xlabel("Epochs", fontsize = 16)
ax.set_ylabel("Cos Similarity", fontsize = 16)

# Set the title
ax.set_title("Figure for Epochs and Cosine Similarities", fontsize=18,backgroundcolor='#F8CECC',
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