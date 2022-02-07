#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 14:29:36 2022

@author: zwl0331
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.io

#%% Xray
Xray = pd.read_pickle(r'Spectrums/XraySpectra/20180409r003s003.pickle')

X_energy = Xray['E'][0]
Xdistribution = Xray['normalised_number']
Xuncertain = Xray['normalised_number_sem']
bin_widths = Xray['bin_widths']

plt.plot(X_energy[:457], Xdistribution[:457])
plt.xlabel('X-ray energy/eV')
plt.ylabel('Normalised number')

#%% Gammaray
gamma = scipy.io.loadmat('Spectrums/GammaSpectra/DatasetA/DatasetA_GammaSpecFits.mat')

g_energy = gamma['GammaEnergy_MeV'][0]
g_number = gamma['ExpNph'][0]
plt.plot(g_energy, g_number)
plt.yscale('log')

