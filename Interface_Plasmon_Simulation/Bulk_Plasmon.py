# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 19:08:06 2017

@author: Eric Hoglund
"""

from __future__ import division, print_function
import numpy as np

hbar=6.582E-16#[eV s]
e=1.602E-19 #[C]
a0 = 5.292*10**-11 #[m]
c=3E8#[m/s]
m0=9.11E-31#[kg]
m0c2=511.#[keV]
eps0=8.85E-12#[C^2/N m]


def bulk_plasmon_double_differential_cross_section(microscope,spectrum,material):
    
    material.set_Ep(q=spectrum.q_perp)
    material.set_eps(E=spectrum.E)
    
    A = 1/(np.pi * a0 * m0c2 * (microscope.v/c)**2 * material.na) #[m^2/eV]
    fdE = np.imag(-1/material.eps)
    B = 1/((spectrum.q_perp/microscope.k0)**2 + (spectrum.E/(hbar*microscope.v)))**2
    return A*fdE*B