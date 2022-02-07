# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 13:17:02 2022

@author: 2dong
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from function import *

#%% Xray
Xray = pd.read_pickle(r'XrayBath/XraySpectra/20180409r003s003.pickle')

X_energy = Xray['E'][0]
Xdistribution = Xray['normalised_number']
Xuncertain = Xray['normalised_number_sem']
bin_widths = Xray['bin_widths']

#%% Gammaray
gamma = scipy.io.loadmat('GammaSpectra/DatasetA/DatasetA_GammaSpecFits.mat')

g_energy = gamma['GammaEnergy_MeV'][0]
g_number = gamma['ExpNph'][0]

#%%
#finding rate
#x-ray radius
rx = c*tx
#interaction angle range
tmin = np.arcsin(d/rx)
tmax = np.pi - tmin
# range of interaction angles
theta = np.linspace(tmin, tmax, 500)
theta_bw = np.array([(tmax - tmin)/499]*500)

#get nx
#filter in range 1.3-1.5 keV
lower_x = np.digitize(1300, X_energy)
upper_x = np.digitize(1500, X_energy, True)
#number of xrays
Nx = Xdistribution[lower_x: upper_x]
#energy of xrays
Ex = X_energy[lower_x: upper_x]
#number density
n = x_density(Nx, theta)

#get Ng
#filter in range 660-760 MeV
lower_g = np.digitize(660, g_energy)
upper_g = np.digitize(760, g_energy, right = True)
#number of gamma rays
Ng = g_number[lower_g: upper_g]
#energy of gamma rays
Eg = g_energy[lower_g: upper_g]

print(find_pairs(theta, theta_bw, Ex, Eg, n, Ng))



