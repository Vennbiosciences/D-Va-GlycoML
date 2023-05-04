#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read epochs, cos_similarities data, and plot them.
The 15 training files are from Energy-01-HCD-15-20-34-37-40 .
The 4 testing file are from
    Energy-01-HCD-15-20-34-37-40  Energy-02-HCD-15-20-35-40  Energy-03-HCD-20-25-40-45  Energy-04-HCD-30.
The batch_size = 8.

Created on 17 April 2022.
###################################################################################################################
"""
__author__ = 'ZLiang'

import matplotlib.pyplot as plt
import numpy as np

# Cosine similarities for different test cases. Totally 34 scores for each test.
# These are used for the y-axis.
test_1_mean = [0.924370, 0.934690, 0.941953, 0.945356, 0.946473, 0.947832, 0.949238, 0.949706, 0.950814, 0.951501,
               0.951591, 0.952122, 0.952589, 0.952236, 0.952632, 0.952517, 0.953037, 0.953628, 0.953245, 0.953362,
               0.953520, 0.954029, 0.953786, 0.953944, 0.953531, 0.953898, 0.954066, 0.954378, 0.954429, 0.954155,
               0.953974, 0.954168, 0.953711, 0.954282, ]
test_1_median = [0.963244, 0.970500, 0.976539, 0.978251, 0.979320, 0.980713, 0.979787, 0.980811, 0.982511, 0.982191,
                 0.982257, 0.982328, 0.982549, 0.982819, 0.983529, 0.982722, 0.984188, 0.983632, 0.983864, 0.983869,
                 0.983656, 0.983801, 0.983186, 0.984468, 0.982545, 0.983873, 0.984597, 0.984414, 0.982611, 0.984495,
                 0.984998, 0.984456, 0.984699, 0.985530, ]
test_2_mean = [0.916296, 0.920391, 0.925850, 0.923536, 0.925096, 0.923131, 0.919353, 0.920033, 0.920205, 0.919212,
               0.919508, 0.921014, 0.918548, 0.918881, 0.920412, 0.919157, 0.918926, 0.918096, 0.917989, 0.915303,
               0.918665, 0.917418, 0.916582, 0.917328, 0.915943, 0.916414, 0.915440, 0.917037, 0.913307, 0.918266,
               0.916093, 0.915663, 0.915698, 0.916308, ]
test_2_median = [0.952597, 0.953405, 0.961162, 0.960534, 0.963109, 0.959994, 0.956462, 0.958313, 0.959546, 0.959148,
                 0.958793, 0.958962, 0.957625, 0.958624, 0.959850, 0.958432, 0.959590, 0.958916, 0.957923, 0.955542,
                 0.959702, 0.958206, 0.957391, 0.958471, 0.956083, 0.956095, 0.956059, 0.957427, 0.953515, 0.958852,
                 0.955120, 0.956418, 0.955739, 0.956738, ]
test_3_mean = [0.844127, 0.830710, 0.834758, 0.839926, 0.831429, 0.838333, 0.837059, 0.838739, 0.831998, 0.832838,
               0.833100, 0.830811, 0.833268, 0.829934, 0.827189, 0.828356, 0.827995, 0.828758, 0.825951, 0.830698,
               0.825350, 0.828954, 0.824814, 0.826886, 0.829308, 0.829197, 0.829583, 0.825714, 0.831927, 0.826316,
               0.829033, 0.826810, 0.829426, 0.827890, ]
test_3_median = [0.887312, 0.869003, 0.876862, 0.882871, 0.876925, 0.876736, 0.875110, 0.880652, 0.875291, 0.875246,
                 0.874710, 0.874701, 0.876358, 0.875678, 0.869912, 0.872724, 0.871008, 0.872478, 0.867599, 0.873616,
                 0.868233, 0.873841, 0.867330, 0.870086, 0.870724, 0.871309, 0.874060, 0.867985, 0.873005, 0.870576,
                 0.868008, 0.869344, 0.872784, 0.870458, ]
test_4_mean = [0.583467, 0.566392, 0.545388, 0.553989, 0.543296, 0.563209, 0.558880, 0.559514, 0.537709, 0.542823,
               0.550117, 0.543475, 0.547815, 0.543936, 0.536092, 0.544323, 0.537579, 0.548099, 0.539944, 0.543692,
               0.543272, 0.531273, 0.536144, 0.532051, 0.547182, 0.544994, 0.542813, 0.535207, 0.547161, 0.536251,
               0.533977, 0.543800, 0.535249, 0.532679, ]
test_4_median = [0.615004, 0.588633, 0.572332, 0.574591, 0.562994, 0.581349, 0.579466, 0.575123, 0.542870, 0.547395,
                 0.564480, 0.554096, 0.562273, 0.559072, 0.537283, 0.552477, 0.538375, 0.560844, 0.553032, 0.546202,
                 0.548986, 0.524430, 0.538557, 0.533312, 0.547960, 0.545980, 0.549779, 0.537153, 0.545009, 0.539225,
                 0.533389, 0.542781, 0.536680, 0.537395, ]


# The number of folds for the training, each fold contains 5 epochs.
# The folds range from 1 to 35, corresponding the epochs range from 5 to 500.
x = np.arange(1, 35)

# Create the size of the figure.
fig=plt.figure(num=1,figsize=(20, 9))

# Set the size of subplot and starting point, set three styles for the plot
ax=fig.add_subplot(111)
ax.plot(x, test_1_mean, "r--.", label="Testing the Mean for Energy-01-HCD-15-20-34-37-40")
ax.plot(x, test_1_median, "r-x", label="Testing the Median for Energy-01-HCD-15-20-34-37-40")
ax.plot(x, test_2_mean, "g-.^", label="Testing the Mean for Energy-02-HCD-15-20-35-40")
ax.plot(x, test_2_median, "g:+", label="Testing the Median for Energy-02-HCD-15-20-35-40")
ax.plot(x, test_3_mean, "b--+", label="Testing the Mean for Energy-03-HCD-20-25-40-45")
ax.plot(x, test_3_median, "b-d", label="Testing the Median for Energy-03-HCD-20-25-40-45")
ax.plot(x, test_4_mean, "y-.s", label="Testing the Mean for Energy-04-HCD-30")
ax.plot(x, test_4_median, "y:p", label="Testing the Median for Energy-04-HCD-30")

# Set the scales for x-axis and y-axis.
ax.set_xlim([1, 34.1])
ax.set_ylim([0, 1])

# Set the scales for the display of x-axis and y-axis.
# np.linspace() is arithmetic progression， from 1 to 34
ax.set_xticks(np.linspace(1, 34, 34))
ax.set_yticks(np.linspace(0, 1, 11))

# Set the labels for the scales
ax.set_xticklabels(["5", "20", "35", "50", "65", "80", "95", "110", "125", "140", "155", "170", "185", "200", "215",
                    "230", "245", "260", "275", "290", "305", "320", "335", "350", "365", "380", "395", "410", "425",
                    "440", "455", "470", "485", "500", ], fontproperties="monospace", fontsize=12, rotation=10)
ax.set_yticklabels(["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])

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