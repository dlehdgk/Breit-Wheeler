# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 13:58:01 2022

@author: 2dong
"""

# need a graph of cross section data for given E_g and E_x
# scalar field display with axes E_g and E_x
# 1. function for sigma (total cross-section of BW for all solid angles) taking E_g, E_x and theta

#imports
import numpy as np
import matplotlib.pyplot as plt

#constants
# electron rest mass in eV
me = 0.511e6
# classical electron radius in m
re = 2.82e-15
#%% function for cross-section

def sigma(E_g, E_x, c):
    # CoM invariant mass^2
    s = 2*(1-c)*E_g*E_x
    if s < 4*me*me:
        sig = 0
    else:
        #beta
        b = np.sqrt(1-(4*me*me)/s)
        sig = 0.5*np.pi*re*re*(1-b*b)*((3-b*b*b*b)*np.log((1+b)/(1-b))-2*b*(2-b*b))
    return sig
#%% datasets for E_g and E_x
# 2. set of data for E_g and E_x with reasonable range in 1D array form

# x-ray spectrum of range 1.3-1.4 keV with resolution 4 eV
E_x = np.linspace(1.3e3, 1.4e3, 100)
# gamma ray maximum of 710 +- 50 MeV
E_g = np.linspace(660e6, 760e6, 100)
#%%
# 3. matrix of sigma values for E_g and E_x values
#matrix of cross-sections 
sig = np.zeros((len(E_g), len(E_x)))

#angle of interaction
cos = -0.1

for i in range(len(E_g)):
    for j in range(len(E_x)):
        sig[i][j] = sigma(E_g[i], E_x[j], cos)
#%% print matrix
#scaling
vmax = sig.max()
vmin = sig.min()
xmin = E_x.min()
xmax = E_x.max()
gmin = E_g.min()
gmax = E_g.max()

x,y = np.meshgrid(E_x, E_g)

plt.title("cross-section for given photon E at cos(theta) = %f" % cos, y=1.08)
plt.xlabel("X-ray energy / eV")
plt.ylabel("Gamma ray energy / eV")

contours = plt.contour(x, y, sig, 6, vmax= vmax, vmin= vmin, colors='black')
plt.clabel(contours, inline=True, fontsize=8)
plt.imshow(sig, extent=[xmin, xmax, gmin, gmax], origin='lower',
           cmap='jet', aspect='auto')
plt.colorbar()  
