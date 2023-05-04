#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, loss entropy data, and plot them.
The 285 training files are from Training-01-Human-285 .
The 285 + 76 training from Training-01-Human-285 and Training-02-Cyno-76 .
The 285 + 155 training from Training-01-Human-285 and Training-03-Human-for-evaluation-155 .

Created on 17 May 2022.
###################################################################################################################
"""
__author__ = 'ZLiang'

import matplotlib.pyplot as plt
import numpy as np

# Cosine similarities for different test cases. Totally 34 scores for each test.
# These are used for the y-axis.
# The values are from the folds of
#   0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48,
#   51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99
train_1_loss_entropy=[3.03869366645813, 2.952838659286499, 2.9502334594726562, 2.9489660263061523, 2.948162794113159,
                      2.9475839138031006, 2.9471077919006348, 2.9467904567718506, 2.9464995861053467, 2.946242570877075,
                      2.9460272789001465, 2.9458210468292236, 2.945687770843506, 2.945544958114624, 2.9454119205474854,
                      2.945291757583618, 2.9452016353607178, 2.945114850997925, 2.945013999938965, 2.9449057579040527,
                      2.944851875305176, 2.9447784423828125, 2.944695234298706, 2.9446380138397217, 2.9445741176605225,
                      2.9445483684539795, 2.9444875717163086, 2.944472074508667, 2.9444570541381836, 2.9443840980529785,
                      2.9443140029907227, 2.9443013668060303, 2.944230318069458, 2.944199562072754]
test_1_median = [0.9610788775558861, 0.9594953760033962, 0.958979173033334, 0.9577877678634839, 0.960146661146974,
                 0.9591893947395627, 0.9573678515992754, 0.9580576507721082, 0.9572387404941064, 0.9584127025751459,
                 0.9582270626831269, 0.9565091841702549, 0.9567606154342047, 0.9575744818355216, 0.9575036183558338,
                 0.9572184627637984, 0.9567079589823593, 0.9568761222682003, 0.957032034618683, 0.9559853101077711,
                 0.9573807116342606, 0.9570235699692926, 0.9567952535829239, 0.959079869663946, 0.9596032852244832,
                 0.9565126566598032, 0.9581063977151558, 0.958148931656986, 0.9570519155981115, 0.9579779553339817,
                 0.9576534537127527, 0.9557151790794052, 0.9564103321479606, 0.9581945368739618, ]

# The number of folds for the training, each fold contains 5 epochs.
# The folds range from 1 to 35, corresponding the epochs range from 5 to 500.
x = np.arange(1, 35)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(16, 8))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, train_1_loss_entropy, "r--.", label="Loss entropy on Training-01")
ax.plot(x, test_1_median, "b-.^", label="Median for Energy-01-HCD-15-20-34-37-40 on Training-01")

# Set the scales for x-axis and y-axis.
ax.set_xlim([1, 34.1])
ax.set_ylim([0.5, 3.1])

# Set the scales for the display of x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 34
ax.set_xticks(np.linspace(1, 34, 34))
ax.set_yticks(np.linspace(0.5, 3.1, 6))

# Set the labels for the scales
ax.set_xticklabels(["5", "20", "35", "50", "65", "80", "95", "110", "125", "140", "155", "170", "185", "200", "215",
                    "230", "245", "260", "275", "290", "305", "320", "335", "350", "365", "380", "395", "410", "425",
                    "440", "455", "470", "485", "500", ], fontproperties="monospace", fontsize=12, rotation=10)
ax.set_yticklabels(["0.5", "1.0", "1.5", "2.0", "2.5", "3.0"])

# Set the decorations for the scales and labels
ax.tick_params(left=False,pad=8,direction="in",length=2,width=3,color="b",labelsize=12)
ax.tick_params("x",labelrotation=10)

# Set the labels for x and y
ax.set_xlabel("Epochs", fontsize = 16)
ax.set_ylabel("Cos Similarity or Loss Entropy", fontsize = 16)

# Set the title
ax.set_title("Figure for Epochs, Cosine Similarities and Loss Entropies", fontsize=18,backgroundcolor='#F8CECC',
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