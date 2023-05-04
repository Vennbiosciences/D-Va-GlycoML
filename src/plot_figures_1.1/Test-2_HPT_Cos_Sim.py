#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###################################################################################################################
This script read a csv file, which contains "cos_all", "cos_by", and "cos_Y". Then plot them separately.

The experimental spectra are from Test-Human-N-High-SCE, and the predicted spectra are from HPT.

Created on 10 January 2023.
###################################################################################################################
"""
__author__ = 'ZLiang'

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# setting
pd.options.display.notebook_repr_html=False  # notebook display
plt.rcParams['figure.dpi'] = 300  # figure dpi
sns.set_theme(style='darkgrid')  # figure style

cos_sim = pd.read_csv(r'Spe-Match-cos-sim-Test-02-HPT.txt', sep='\t', encoding='utf-8')
cos_sim.head()

# x-axis histogram for 'Cos_All_0.001'
#sns.histplot(data=cos_sim, x='Cos_All_HPT_Test-Human-N-High-SCE', bins=100, color='b', kde=True)

# x-axis histogram for 'Cos_by_0.001'
#sns.histplot(data=cos_sim, x='Cos_by_HPT_Test-Human-N-High-SCE', bins=100, color='r', kde=True)

# x-axis histogram for 'Cos_All_0.001'
sns.histplot(data=cos_sim, x='Cos_Y_HPT_Test-Human-N-High-SCE', bins=100, color='g', kde=True)

plt.show()