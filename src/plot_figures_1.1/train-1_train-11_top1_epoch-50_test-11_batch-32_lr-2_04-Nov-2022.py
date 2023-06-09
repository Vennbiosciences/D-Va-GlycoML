#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The training files Training-01-all21-N-Human345_Cyno76-Full_Semi.
And use transfer learning for Training-11-all21-O-Human345

The testing files are from Testing-11-O-2202-P4M .

The batch size is 32, learning rates are 2.

Created on 04 November 2022.
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


test_11_train_11_mean=[0.8769256418855105, 0.8836770114749661, 0.8952600492955058, 0.9001049657113833, 0.9115907221130886,
                       0.9035566055514762, 0.9110414822021833, 0.915485103031655,  0.9009513772150919, 0.9132407457649958,
                       0.9091391664520184, 0.9119437934396472, 0.9218883255886782, 0.9138293486767693, 0.9216966790564304,
                       0.9137427273052745, 0.9207493580686607, 0.9235329613178472, 0.9132341048563626, 0.9193821483647442,
                       0.9244014029559469, 0.9218269551579458, 0.916478271626511,  0.9182223224266219, 0.9215779200062199,
                       0.9246564342316262, 0.9225735766017977, 0.9263579069110367, 0.91892827948015,   0.9213161531765107,
                       0.9187883883516631, 0.9200124767527728, 0.9187319783441248, 0.9266060321229763, 0.9246330345635424,
                       0.9195897793286445, 0.9215765041719706, 0.9262117276598416, 0.9241367260621899, 0.9240629194589998,
                       0.9262308335081416, 0.9284822758956695, 0.9235315862386632, 0.9241014402695933, 0.9250490288561359,
                       0.9244120587288917, 0.9267108870607864, 0.921801268994611,  0.923210624759032,  0.9271579790495749,
                       ]

test_11_train_11_median=[0.9159432452199527, 0.9211807862784335, 0.9418595603006842, 0.9502499908015405, 0.9472461979845646,
                         0.9378717232850912, 0.9466947583442259, 0.9574451764114335, 0.9568629533961299, 0.9560595495142872,
                         0.958563460190681,  0.9578678936244711, 0.9586874063261228, 0.9557782507130109, 0.961759580487161,
                         0.9555687943417781, 0.9664216681681879, 0.9580219684998186, 0.9610576982629095, 0.9649013919828994,
                         0.9642308803415611, 0.9680344856984406, 0.9683585285567584, 0.957841557224112,  0.9644901137167402,
                         0.9660851745759566, 0.965168852808681,  0.9635731036298552, 0.9634188963677061, 0.9678487875890861,
                         0.9595512206285048, 0.9672093738129834, 0.9614222351279239, 0.9678484129050551, 0.9644510128309227,
                         0.96902183409419,   0.9705352157974023, 0.9693067361965985, 0.9683834911384779, 0.9629488419572174,
                         0.9678347378129051, 0.9666398686509619, 0.9700208593722942, 0.9671182960554883, 0.9670848097124667,
                         0.9705392949260248, 0.9656006814194307, 0.9704518829724305, 0.9695833441111374, 0.9660181748964729,
                         ]

test_11_train_01_11_mean=[0.914773319535638,  0.9214363438179931, 0.9218871675839119, 0.921261396823338,  0.9239595848301593,
                          0.9214164545981404, 0.9257871909729085, 0.9232683061911595, 0.9254568116700332, 0.9242305861052482,
                          0.9200771723012874, 0.924300877918673,  0.9258847988251442, 0.924493105377621,  0.922833982300992,
                          0.9234956286560998, 0.9246542850891544, 0.9244021753600888, 0.9263132930709022, 0.9269948534535629,
                          0.923556304573571,  0.9250232758878726, 0.9242407277588868, 0.9274344504987779, 0.9275160443138895,
                          0.9238065759742675, 0.9245817947553912, 0.9260458951327021, 0.9244108590308131, 0.9252936888918296,
                          0.9241100006332047, 0.9241777816900422, 0.9238697762313416, 0.9232110912127907, 0.9282806211047074,
                          0.9267915484332574, 0.9244954464954835, 0.9251419495670524, 0.9255225014168338, 0.9269158940874462,
                          0.9262433312905992, 0.9238614307829982, 0.9230718068578879, 0.9236604996972647, 0.9265362534322241,
                          0.9240714340138322, 0.9250924377649784, 0.9281215835090301, 0.9261286772426792, 0.9249895604608929,
                          ]

test_11_train_01_11_median=[0.9630004585531933, 0.9665058756051046, 0.9663834997756433, 0.9629010322793494, 0.967028317337196,
                            0.9668499415833018, 0.9668120697630538, 0.9652281224268322, 0.9657149566341298, 0.9690953075993645,
                            0.9666139279094872, 0.9687189026753633, 0.9684730454982298, 0.9698437248790153, 0.9697008436711472,
                            0.9694735703658248, 0.9694094253928729, 0.9659760551260823, 0.9695337749560234, 0.9693899338903971,
                            0.9699760562474593, 0.9701453910969018, 0.9700519853655967, 0.9710307667945145, 0.9681330346512704,
                            0.9722389665942859, 0.970438647724489,  0.9694941341856611, 0.9696822879623752, 0.969154065745129,
                            0.9704589650767503, 0.969250053744986,  0.9687690737084382, 0.9713893193202553, 0.9684323993349238,
                            0.9701539665051637, 0.9699084278414027, 0.9709862802102698, 0.9698030242575806, 0.9677978654738271,
                            0.9690293095665625, 0.9699369878156907, 0.9706031080106838, 0.970660936588164,  0.9693537750278158,
                            0.969823525948273,  0.9688946759897152, 0.9673260088422045, 0.9692534022841894, 0.9706924643453382,
                            ]

# The number of epochs for the training.
x = np.arange(1, 51)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(20, 10))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_11_train_11_mean, "y--.", label="Mean for Testing-11 on Training-11, batch size=32, learning rate=0.0002")
ax.plot(x, test_11_train_11_median, "g:+", label="Median for Testing-11 on Training-11, batch size=32, learning rate=0.0002")

ax.plot(x, test_11_train_01_11_mean, "b--+", label="Mean for Testing-11 on transfer learning Training-01, batch size=32, learning rate=0.0002")
ax.plot(x, test_11_train_01_11_median, "r-x", label="Median for Testing-11 on transfer learning Training-01, batch size=32, learning rate=0.0002")

# Set the scales for x-axis and y-axis.
ax.set_xlim([1, 50.1])
ax.set_ylim([0.86, 0.98])

# Set the scales for the display of x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 50
ax.set_xticks(np.linspace(1, 50, 50))
ax.set_yticks(np.linspace(0.86, 0.98, 4))

# Set the labels for the scales
ax.set_xticklabels(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
                    "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34",
                    "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50"],
                   fontproperties="monospace", fontsize=12, rotation=10)
ax.set_yticklabels(["0.86", "0.90", "0.94", "0.98"])

# Set the decorations for the scales and labels
ax.tick_params(left=False,pad=8,direction="in",length=2,width=3,color="b",labelsize=12)
ax.tick_params("x",labelrotation=10)

# Set the labels for x and y
ax.set_xlabel("Epochs", fontsize = 16)
ax.set_ylabel("Cos Similarity", fontsize = 16)

# Set the title
ax.set_title("Figure for Epochs and Cosine Similarities, Training-01, Training-11, Testing-11", fontsize=18,backgroundcolor='#DAE8FC',
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