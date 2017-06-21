# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 17:52:20 2017

@author: Eric Hoglund
"""

from __future__ import division, print_function
from traits.api import HasTraits, Float, List, Instance, on_trait_change
from traitsui.api import *
import numpy as np

hbar=6.582E-16#[eV s]
e=1.602E-19 #[C]
a0 = 5.292*10**-11 #[m]
c=3E8#[m/s]
m0=9.11E-31#[kg]
m0c2=511.#[keV]
eps0=8.85E-12#[C^2/N m]

class microscope(HasTraits):
    #__parent = Instance(system)
    def __init__(self, keV=100, resolution=0.05, collection_angle=2E-3):
        parent_system = None
        self.keV = keV #self.keV={'value':keV, 'units':'keV'}
        
        self.gamma = 1 + keV/m0c2 #[au]
        
        v2_c2 = self.keV * (keV + 2 * m0c2)/(keV + m0c2)**2
        self.v = c * np.sqrt(v2_c2) #[m/s]
        self.T = m0*self.v**2/(2*e)#m0c2/2 * (self.v/c)**2 #[keV]
        self.resolution = resolution #Nion=0.05, Titan=0.05
        self.collection_angle = collection_angle #[mrad]
        
        self.k0 = m0*self.v*self.gamma/(hbar*e) #[1/m]
        
        @on_trait_change('parent_system')
        def _parent_system_changed():
            print('Parent sys',parent_system)
            parent_system.test = 4.0
            if parent_system is not None:
                print('Printing parent_system test',parent_system.test)
    def print_parameters(self):
        print('Microscope')
        for each in self.__dict__:
            print(each,': ',self.__dict__[each])
        print()