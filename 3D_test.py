# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 14:08:44 2022

@author: 2dong
"""

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

# x-ray spectrum of range 1.3-1.4 keV with resolution 4 eV y axis
E_x = np.linspace(1e3, 200e3, 200)
# gamma ray maximum of 710 +- 50 MeV z axis
E_g = np.linspace(110e6, 1000e6, 100)
cos = 0.3

#%%
# 3. matrix of sigma values for E_g and E_x values
#matrix of cross-sections 
sig = np.zeros((len(E_g), len(E_x)))

for i in range(len(E_g)):
    for j in range(len(E_x)):
        sig[i][j] = sigma(E_g[i], E_x[j], cos)
x,y = np.meshgrid(E_x, E_g)
#no.gamma per shot
N_g = 7e7
#number density of x-ray
n_x = 1.4e21
# diameter of x-ray cloud
z = 97e-12 * 3e8
#rate of BW in thin target approx with const no. gamma and number density of x-ray
R = N_g * n_x * z * sig
fig = plt.figure()
ax = plt.axes(projection='3d')
p = ax.plot_surface(x, y, R, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('rough BW rate for cos(theta) = %f' % cos)
ax.set_xlabel('x-ray')
ax.set_ylabel('$\gamma$ ray')
ax.set_zlabel('rate per shot')
fig.colorbar(p, pad = 0.1)











