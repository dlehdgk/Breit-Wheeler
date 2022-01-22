# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 15:59:53 2022

@author: 2dong

checking which angles of theta allow for the BW for given energies
"""

#imports
import numpy as np
import matplotlib.pyplot as plt

#constants
# electron rest mass in eV
me = 0.511e6
#%% in cos theta
def s(E, c):
    s = 2*(1-c)*E* 710e6
    if s < 4*me*me:
        s = 0
    return s
#x-ray spread for fixed gamma
E = np.linspace(1.3e3, 1.4e3, 100)
c = np.linspace(-1, 1, 100)

m = np.zeros((len(E), len(c)))

for i in range(len(E)):
    for j in range(len(c)):
        m[i][j] = s(E[i], c[j])
        
x,y = np.meshgrid(c, E)
#%% graph scales
ma = m[m != 0]
vmax = ma.max()
vmin = ma.min()
Emin = E.min()
Emax = E.max()
#%%
plt.title("invariant mass^2 for x-ray and angle with gamma = 710 MeV", y=1.08)
plt.xlabel("cos(interaction angle)")
plt.ylabel("x-ray energies / eV")
contours = plt.contour(x, y, m, 3, vmax= vmax, vmin= vmin, colors='black')
plt.clabel(contours, inline=True, fontsize=8)
plt.imshow(m, extent=[-1, 1, Emin, Emax], origin='lower',
           cmap='jet', aspect='auto')
plt.colorbar()