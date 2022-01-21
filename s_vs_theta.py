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

#E = E_g x E_x
def s(E, t):
    s = 2*(1-np.cos(t))*E
    if s < 4*me*me:
        s = 0
    return s
#x-ray spread for fixed gamma
E = np.linspace(1.3e3, 1.4e3, 50) * 710e6
theta = np.linspace(0, np.pi, 50)

m = np.zeros((len(E), len(theta)))

for i in range(len(E)):
    for j in range(len(theta)):
        m[i][j] = s(E[i], theta[j])
        
x,y = np.meshgrid(theta, E)

plt.title("s for given E^2 and angle")
plt.xlabel("interaction angle / rads")
plt.ylabel("710MeV x Exray / eV^2")
plt.contourf(x, y, m, cmap = 'jet')
plt.colorbar()
plt.show()
