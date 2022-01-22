# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 09:26:54 2022

@author: 2dong
"""
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
#%%
#plotting cross-section for range 0<theta<pi and 1.3 keV < x-ray < 1.4 keV and 710 MeV gamma
#E gamma = 710 MeV
Eg = 710e6
c = np.linspace(-1, 1, 100)
E_x = np.linspace(1.3e3, 1.4e3, 100)
# x-section for Ex and theta
sig = np.zeros((len(E_x), len(c)))

for i in range(len(E_x)):
    for j in range(len(c)):
        sig[i][j] = sigma(Eg, E_x[i], c[j])

x,y = np.meshgrid(c, E_x)

ma = sig[sig != 0]
vmax = ma.max()
vmin = ma.min()

Emin = E_x.min()
Emax = E_x.max()
#%%
plt.title("cross-section for given x-ray energies and angles", y=1.08)
plt.xlabel("cos(interaction angle)")
plt.ylabel("x-ray energy / eV")
contours = plt.contour(x, y, sig, 4, vmax= vmax, vmin= vmin, colors='black')
plt.clabel(contours, inline=True, fontsize=8)

plt.imshow(sig, extent=[-1, 1, Emin, Emax], origin='lower',
           cmap='jet', aspect='auto')
plt.colorbar()
        
