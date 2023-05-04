#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The three training files are from Test 1, Test 2, and Test 3.
The testing file are from Test 3.
The batch_size = 8 .

Created on 29 March 2022.
###################################################################################################################
"""
__author__ = 'ZLiang'

import matplotlib.pyplot as plt
import numpy as np

# Cosine similarities for different test cases. Totally 19 scores for each test.
# These are used for the y-axis.
test_1 = [0.421085, 0.461061, 0.469588, 0.443070, 0.440987, 0.465916, 0.474323, 0.452174, 0.454342, 0.453191, 0.479363,
          0.473990, 0.452135, 0.468929, 0.460150, 0.466591, 0.455720, 0.445931, 0.466425, ]
test_2 = [0.630641, 0.599521, 0.617129, 0.625572, 0.620913, 0.627726, 0.626242, 0.642002, 0.626551, 0.599636, 0.624680,
          0.609387, 0.606494, 0.601542, 0.602350, 0.622688, 0.613510, 0.621278, 0.614985, ]
test_3 = [0.922272, 0.933559, 0.943946, 0.934209, 0.944586, 0.950660, 0.951656, 0.953720, 0.954094, 0.957394, 0.956975,
          0.958325, 0.958726, 0.958988, 0.957584, 0.959192, 0.960877, 0.961285, 0.961804, ]

# The number of folds for the training, each fold contains 5 epochs.
# The folds range from 1 to 20, corresponding the epochs range from 5 to 95.
x = np.arange(1, 20)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(15,8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1, "g-.d", label="Training on Test 1, Testing on Test 3")
ax.plot(x, test_2, "r-x", label="Training on Test 2, Testing on Test 3")
ax.plot(x, test_3, "b:+", label="Training on Test 3, Testing on Test 3")

# Set the scales for x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 19
ax.set_xlim([1, 19.1])
ax.set_ylim([0, 1])

# Set the scales for the display of x-axis and y-axis.
ax.set_xticks(np.linspace(1, 19, 19))
ax.set_yticks(np.linspace(0, 1, 11))

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
ax.set_title("Figure for Epochs and Cosine Similarities", fontsize=18,backgroundcolor='#DAE8FC',
             fontweight='bold',color='black',verticalalignment="baseline")

# Set the spines
ax.spines["left"].set_color("darkblue")
ax.spines["bottom"].set_linewidth(2)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

# Set the annotations
ax.annotate(s="Max:0.961804",xy=(19, 0.94),xytext=(18, 0.87),arrowprops=dict(facecolor="b",shrink=0.05,headwidth=12,
                                                                          headlength=6, width=4),fontsize=12)
ax.annotate(s="Min:0.421085",xy=(1,0.4),xytext=(1.2,0.33),arrowprops=dict(facecolor="g",shrink=0.05,headwidth=12,
                                                                          headlength=6, width=4),fontsize=12)

# Set the legend
ax.legend(loc=3,labelspacing=1,handlelength=3,fontsize=14,shadow=True)

plt.show()