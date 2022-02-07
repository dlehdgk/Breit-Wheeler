# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 15:21:33 2022

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
#laser energy (p6 BW paper)
E_l = 2e15 * 40e-12 /100
c = 3e8
#duration of x-ray pulse (p7 BW paper)
tx = 97e-12
#perpendicular interaction distance from Germanium target to gamma trajectory. The default is 1e-3 m (p4 BW paper).
d = 1e-3
#%%
#get list of N_x for given E_x
#1) N_x * laser energy (5J) * 4pi = total number of x-ray photons released in 40 ps pulse
#2) convert to number density for given angle
# - need binned angles from theta min ~ max 
# - nested for loop to run over all angles and energies to get matrix of n_x

def x_density(N_x, theta):
    """
    Parameters
    ----------
    N_x : Number of x-ray photons / eV
        array of the number of x-ray photons per J of laser per srad
    theta : interaction angle / rads
        array of interaction angles of x-ray field and gamma ray pulse

    Returns
    -------
    n : number density matrix
        number density of x-ray field [x-ray energy index, angle index]
    """
    #total number of x-ray photons
    N = N_x * E_l
    #matrix of n[energies, angles]
    n = np.zeros((len(N), len(theta)))
    for i in range (len(N)):
        for j in range (len(theta)):
            n[i, j] = N[i]*np.sin(theta[j])*np.sin(theta[j])/(d*d*c*tx)
    return n
#%%
#cross-section
def sigma(E_g, E_x, theta):
    c = np.cos(theta)
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
#BW pairs per shot
def find_pairs(theta, t_bw, Ex, Eg, nx, Ng):
    """
    Parameters
    ----------
    theta : interaction angle / rads
    
    t_bw : theta bandwidth
    
    Ex : x-ray energies / eV
    
    Ex_bw : x-ray energy bandwidth
    
    Eg : $\gamma$ ray energies / MeV
    
    Eg_bw : $\gamma$ ray energy bandwidth
    
    nx : x-ray number density [x-ray index, angle index] / m-3

    Ng : number of $\gamma$ ray photons per shot

    d : interaction distance / m
        perpendicular interaction distance from Germanium target 
        to gamma trajectory. The default is 1e-3 m (p4 BW paper).

    Returns
    -------
    pairs : number of pairs per shot

    """
    pairs = 0
    for i in range(len(theta)):
        for j in range(len(Ex)):
            for k in range(len(Eg)):
                pairs += d*sigma(Eg[k] * 1e6, Ex[j], np.cos(theta[i])) \
                    *Ng[k]*nx[j,i]*t_bw[i]/(np.sin(theta[i])*np.sin(theta[i]))
    return pairs
