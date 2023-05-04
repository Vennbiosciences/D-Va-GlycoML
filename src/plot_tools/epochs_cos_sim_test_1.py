#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The three training files are from Test 1, Test 2, and Test 3.
The testing file are from Test 1.
The batch_size = 8 .

Created on 29 March 2022.
###################################################################################################################
"""
__author__ = 'ZLiang'

import matplotlib.pyplot as plt
import numpy as np

# Cosine similarities for different test cases. Totally 19 scores for each test.
# These are used for the y-axis.
test_1 = [0.896399, 0.915499, 0.927353, 0.931198, 0.937644, 0.942555, 0.944815, 0.947100, 0.948178, 0.949395, 0.952743,
          0.952673, 0.954353, 0.956802, 0.958024, 0.958974, 0.960398, 0.960302, 0.960304, ]
test_2 = [0.764321, 0.804685, 0.799082, 0.795961, 0.813574, 0.808678, 0.810407, 0.800530, 0.818176, 0.810684, 0.816222,
          0.819429, 0.814971, 0.814973, 0.820029, 0.812236, 0.803986, 0.804978, 0.818378, ]
test_3 = [0.466598, 0.423826, 0.436674, 0.428044, 0.427360, 0.435090, 0.422932, 0.429382, 0.436274, 0.432207, 0.436836,
          0.430710, 0.432816, 0.427147, 0.428556, 0.436694, 0.436655, 0.427251, 0.426848, ]

# The number of folds for the training, each fold contains 5 epochs.
# The folds range from 1 to 20, corresponding the epochs range from 5 to 95.
x = np.arange(1, 20)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(15,8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1, "g-.d", label="Training on Test 1, Testing on Test 1")
ax.plot(x, test_2, "r-x", label="Training on Test 2, Testing on Test 1")
ax.plot(x, test_3, "b:+", label="Training on Test 3, Testing on Test 1")

# Set the scales for x-axis and y-axis.
ax.set_xlim([1, 19.1])
ax.set_ylim([0, 1])

# Set the scales for the display of x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 19
ax.set_xticks(np.linspace(1, 19, 19))
ax.set_yticks(np.linspace(0, 1, 11))

# Set the labels for the scales
ax.set_xticklabels(["5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65", "70", "75", "80",
                    "85", "90", "95", ], fontproperties="monospace", fontsize=12, rotation=10)
ax.set_yticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])

# Set the decorations for the scales and labels
ax.tick_params(left=False,pad=8,direction="in",length=2,width=3,color="b",labelsize=12)
ax.tick_params("x",labelrotation=10)

# Set the labels for x and y
ax.set_xlabel("Epochs", fontsize = 16)
ax.set_ylabel("Cos Similarity", fontsize = 16)

# Set the title
ax.set_title("Figure for Epochs and Cosine Similarities", fontsize=18,backgroundcolor='#D5E8D4',
             fontweight='bold',color='black',verticalalignment="baseline")

# Set the spines
ax.spines["left"].set_color("darkblue")
ax.spines["bottom"].set_linewidth(2)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

#ax.text(17,1,"max:0.960398",fontsize=14,color="g",alpha=1)
#ax.text(6,86,"max:90",fontsize=12,alpha=1)

# Set the annotations
ax.annotate(s="Max:0.960398",xy=(17,0.94),xytext=(17.2,0.87),arrowprops=dict(facecolor="g",shrink=0.05,headwidth=12,
                                                                          headlength=6, width=4),fontsize=12)
ax.annotate(s="Min:0.422932",xy=(7,0.4),xytext=(7.2,0.33),arrowprops=dict(facecolor="b",shrink=0.05,headwidth=12,
                                                                          headlength=6, width=4),fontsize=12)

# Set the legend
ax.legend(loc=3, labelspacing=1, handlelength=3, fontsize=14, shadow=True)

plt.show()