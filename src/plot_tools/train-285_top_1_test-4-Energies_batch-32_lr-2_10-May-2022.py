#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 285 training files are from Training-01-Human-285 .
The 4 testing file are from Testing-01-Different-HCD-energies-24:
    (Energy-01-HCD-15-20-34-37-40  Energy-02-HCD-15-20-35-40  Energy-03-HCD-20-25-40-45  Energy-04-HCD-30)
The batch size is 32, learning rate is 2.

Created on 10 May 2022.
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
test_1_median = [0.9610788775558861, 0.9594953760033962, 0.958979173033334, 0.9577877678634839, 0.960146661146974,
                 0.9591893947395627, 0.9573678515992754, 0.9580576507721082, 0.9572387404941064, 0.9584127025751459,
                 0.9582270626831269, 0.9565091841702549, 0.9567606154342047, 0.9575744818355216, 0.9575036183558338,
                 0.9572184627637984, 0.9567079589823593, 0.9568761222682003, 0.957032034618683, 0.9559853101077711,
                 0.9573807116342606, 0.9570235699692926, 0.9567952535829239, 0.959079869663946, 0.9596032852244832,
                 0.9565126566598032, 0.9581063977151558, 0.958148931656986, 0.9570519155981115, 0.9579779553339817,
                 0.9576534537127527, 0.9557151790794052, 0.9564103321479606, 0.9581945368739618, ]
test_2_median = [0.9419699211113807, 0.9428548367354384, 0.9402257487801733, 0.938526775251687, 0.9409378707818649,
                 0.9395589778962129, 0.9375114480535042, 0.9407293239909739, 0.9375996583818222, 0.9390933762626514,
                 0.9406503458106157, 0.9354579403546323, 0.9390651266239539, 0.9377901519418231, 0.9395562419819061,
                 0.9374614306578791, 0.9363963743476603, 0.9388165031510067, 0.9371607992261465, 0.9398075891595046,
                 0.9380662963520072, 0.9367643716663894, 0.9376100107611425, 0.9405052889598358, 0.9415062526069629,
                 0.9363815176491492, 0.9383331116955896, 0.9404258712025387, 0.939185144551397, 0.9381224981006718,
                 0.9379201145972321, 0.9371439765894272, 0.9369873720616768, 0.9367138455887054, ]
test_3_median = [0.8464842280292842, 0.8471510180167758, 0.8416685736917611, 0.8521571728484992, 0.8432955333294303,
                 0.8420982240642654, 0.8459182934550719, 0.8399749415661909, 0.8449674537145568, 0.8406192535836914,
                 0.844618593979986, 0.8430536928672586, 0.8391598373202798, 0.8416062872389922, 0.8381600259573101,
                 0.8407582257697235, 0.8398356742171365, 0.8395419540609917, 0.8442582030494262, 0.8380728774244985,
                 0.8419373396621945, 0.8442214285233118, 0.8338537851465941, 0.8416685736917611, 0.8372755074560139,
                 0.8365525227070564, 0.8391232609565517, 0.8406904212789594, 0.8347728067399347, 0.8411842001063774,
                 0.8402706995963181, 0.8358890651723258, 0.8417674108242247, 0.8311864148017372, ]
test_4_median = [0.5899086081117211, 0.5673761985048669, 0.5807647457437528, 0.5957570241612788, 0.5761681238741134,
                 0.581827829414423, 0.5791529753624438, 0.569463305784446, 0.5749024060150103, 0.5774880371996616,
                 0.555559495869341, 0.5865741286919035, 0.5792173795562496, 0.5764444755839909, 0.5668776992094433,
                 0.571359176203067, 0.5650817972818842, 0.5784808491637812, 0.5679621173925198, 0.5535745616331884,
                 0.5460948189593111, 0.5703159413170535, 0.5640467732742658, 0.5667333604181795, 0.5538326630998399,
                 0.5531441864045685, 0.5543373468189718, 0.5626112575641765, 0.5616765880598482, 0.563532180407079,
                 0.5694782824418427, 0.5644419418667836, 0.567297101860985, 0.5659239065548176, ]


# The number of folds for the training, each fold contains 5 epochs.
# The folds range from 1 to 35, corresponding the epochs range from 5 to 500.
x = np.arange(1, 35)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(20, 9))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1_median, "r-x", label="Testing the Median for Energy-01-HCD-15-20-34-37-40")
ax.plot(x, test_2_median, "g:+", label="Testing the Median for Energy-02-HCD-15-20-35-40")
ax.plot(x, test_3_median, "b-d", label="Testing the Median for Energy-03-HCD-20-25-40-45")
ax.plot(x, test_4_median, "y:p", label="Testing the Median for Energy-04-HCD-30")

# Set the scales for x-axis and y-axis.
ax.set_xlim([1, 34.1])
ax.set_ylim([0.3, 1])

# Set the scales for the display of x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 34
ax.set_xticks(np.linspace(1, 34, 34))
ax.set_yticks(np.linspace(0.3, 1, 8))

# Set the labels for the scales
ax.set_xticklabels(["5", "20", "35", "50", "65", "80", "95", "110", "125", "140", "155", "170", "185", "200", "215",
                    "230", "245", "260", "275", "290", "305", "320", "335", "350", "365", "380", "395", "410", "425",
                    "440", "455", "470", "485", "500", ], fontproperties="monospace", fontsize=12, rotation=10)
ax.set_yticklabels(["0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])

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