# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 19:08:06 2017

@author: Eric Hoglund
"""

from __future__ import division, print_function
import numpy as np
from materials import Al, GB, vac

hbar=6.582E-16#[eV s]
e=1.602E-19 #[C]
a0 = 5.292*10**-11 #[m]
c=3E8#[m/s]
m0=9.11E-31#[kg]
m0c2=511.#[keV]
eps0=8.85E-12#[C^2/N m]

#class microscope():
#    def __init__(self,keV=100, resolution=0.05, collection_angle=2E-3):
#        self.keV=keV
#        self.v=c * np.sqrt(keV * (keV + 2 * m0c2)/(keV + m0c2)**2) #[m/s]
#        self.resolution = resolution #Nion=0.05, Titan=0.05
#        self.collection_angle = collection_angle #[mrad]
#        self.gamma = 1/np.sqrt(1-self.v**2/c**2) #[au]
#        self.T = 1/2*m0c2*(self.v/c)**2 #[keV]
#        self.k0 = m0*self.v*self.gamma/(hbar*e) #[1/m]
#    def print_parameters(self):
#        print('q_beta: %.2f [1/nm]'%(self.k0*self.collection_angle/10**9))
class microscope():
    def __init__(self,keV=100, resolution=0.05, collection_angle=2E-3):
        self.keV=keV
        
        self.gamma = 1 + keV/m0c2 #[au]
        
        self.v=c * np.sqrt(self.keV * (keV + 2 * m0c2)/(keV + m0c2)**2) #[m/s]
        self.T = m0c2/2 * (self.v/c)**2 #[keV]
        self.resolution = resolution #Nion=0.05, Titan=0.05
        self.collection_angle = collection_angle #[mrad]
        
        self.k0 = m0*self.v*self.gamma/(hbar*e) #[1/m]
    def print_parameters(self):
        print('q_beta: %.2f [1/nm]'%(self.k0*self.collection_angle/10**9))

na = 4 * 4.0**3
class bulk_plasmon_double_differential_cross_section():
    def __init__(self,microscope):
        self.A = 1/(np.pi * a0 * m0c2 * (microscope.v/c)**2 * na) #[m^2/eV]
    f_dE = material.set
    
scope = microscope()
test = bulk_plasmon_double_differential_cross_section(microscope())
print(scope.v**2/c**2)
print(scope.gamma)
print(scope.T)
print(scope.k0)

