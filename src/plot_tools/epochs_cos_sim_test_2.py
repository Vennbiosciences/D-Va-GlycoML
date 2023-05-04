#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The three training files are from Test 1, Test 2, and Test 3.
The testing file are from Test 2.
The batch_size = 8 .

Created on 30 March 2022.
###################################################################################################################
"""
__author__ = 'ZLiang'

import matplotlib.pyplot as plt
import numpy as np

# Cosine similarities for different test cases. Totally 19 scores for each test.
# These are used for the y-axis.
test_1 = [0.794387, 0.807395, 0.819946, 0.812569, 0.815159, 0.809013, 0.805569, 0.802476, 0.819866, 0.815612, 0.800674,
          0.810803, 0.805539, 0.816290, 0.795077, 0.823482, 0.807691, 0.804451, 0.816576, ]
test_2 = [0.904688, 0.925091, 0.932445, 0.939708, 0.942204, 0.942614, 0.947577, 0.950003, 0.949890, 0.952793, 0.954191,
          0.955387, 0.954815, 0.957101, 0.959279, 0.960245, 0.959967, 0.961443, 0.962805, ]
test_3 = [0.604427, 0.608106, 0.589960, 0.610871, 0.624513, 0.603911, 0.613802, 0.606974, 0.611726, 0.615085, 0.611343,
          0.611388, 0.606175, 0.606268, 0.616836, 0.619151, 0.608422, 0.617280, 0.612191, ]

# The number of folds for the training, each fold contains 5 epochs.
# The folds range from 1 to 20, corresponding the epochs range from 5 to 95.
x = np.arange(1, 20)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(15,8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1, "g-.d", label="Training on Test 1, Testing on Test 2")
ax.plot(x, test_2, "r-x", label="Training on Test 2, Testing on Test 2")
ax.plot(x, test_3, "b:+", label="Training on Test 3, Testing on Test 2")

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
ax.set_title("Figure for Epochs and Cosine Similarities", fontsize=18, backgroundcolor='#F8CECC',
             fontweight='bold', color='black', verticalalignment="baseline")

# Set the spines
ax.spines["left"].set_color("darkblue")
ax.spines["bottom"].set_linewidth(2)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

#ax.text(7,97,"max:96",fontsize=14,color="g",alpha=1)
#ax.text(6,86,"max:90",fontsize=12,alpha=1)

# Set the annotations
ax.annotate(s="Max:0.962805",xy=(19, 0.94),xytext=(18, 0.87),arrowprops=dict(facecolor="r",shrink=0.05,headwidth=12,
                                                                          headlength=6, width=4),fontsize=12)
ax.annotate(s="Min:0.589960",xy=(3, 0.56),xytext=(3.2,0.49),arrowprops=dict(facecolor="g",shrink=0.05,headwidth=12,
                                                                          headlength=6, width=4),fontsize=12)

# Set the legend
ax.legend(loc=3,labelspacing=1,handlelength=3,fontsize=14,shadow=True)

plt.show()