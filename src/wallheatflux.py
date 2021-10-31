# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 10:18:49 2021

@author: quent
"""

import cantera as ct
import numpy as np

def hFlame(gas,Tw):
    M = gas.mean_molecular_weight
    kb = ct.boltzmann
    Na = ct.avogadro
    N = gas.density_mole
    n = N * Na 
    m = M / Na
    h = n * Tw * (2*kb)**(3/2) / ( 2 * ( np.pi*m*Tw )**(1/2) )
    return(h)

def meanFreePathCalculation(gas):
    mu = gas.viscosity
    M = gas.mean_molecular_weight
    p = gas.P
    T = gas.T
    R = ct.gas_constant
    mfp = mu / p * ( (np.pi*R*T) / (2*M) )**(1/2)
    return(mfp)

gas = ct.Solution('gri30.yaml')

gas.TPX = 300,1e5,'CH4:1,O2:2,N2:{}'.format(2*0.79/0.21)

Twall = 300

h = hFlame(gas,Twall)

mfp = meanFreePathCalculation(gas)

print(h,mfp)