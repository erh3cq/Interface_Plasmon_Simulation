# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 17:52:20 2017

@author: Eric Hoglund
"""

from __future__ import division, print_function
from traits.api import HasTraits, Float, Int, Str, Instance, on_trait_change
from traitsui.api import *
import numpy as np

hbar=6.582E-16#[eV s]
e=1.602E-19 #[C]
a0 = 5.292*10**-11 #[m]
c=3E8#[m/s]
m0=9.11E-31#[kg]
m0c2=511.#[keV]
eps0=8.85E-12#[C^2/N m]

class Microscope(object):
    """
    Microscope class that derives important parameters based on input such as the accelerating voltage.
    
    Adjustable parameters
    ----------------------
    Input parameters that are adjustable.
    
    name : Str
        The name of the microscope.
    keV : Int
        Accelerating voltage.
        Units in keV.
    resolution : Float
        Energy resolution of the microscope (ie FWHM of the ZLP).
        Units in eV.
    dispersion : Float
        Energy dispersion of the spectrometer.
        Units in eV.
    collection_angle : Float
        Spectrometr collection semi-angle.
        Units in mrad.
    
    Derived parameters
    -------------------
    :math:`β^2` : Float
        Relativistic velocity to light.
        :math:`\β^2=v^2/c^2`
    :math:`\gamma` : Float
        Relativistic factor for an electron in vacuum.
        :math:`\gamma^2=1-v^2/c^2`
    v : Float
        Velocity of the electrons in vacuum.
        Units in m/s.
    T : Float
        Kinetic energy of the electrons.
        Units in keV
    :math:`k_0` : Float
        Wave number.
        Units in 1/m.
    """
    
    def __init__(self, keV=100, resolution=0.05, dispersion=0.05, collection_angle=2E-3, **traits):
        name = 'User specified' #desc='Name of the microscope', label='Name')
    
        self.keV = keV #desc='Accelerating voltage', label='keV'
        self.set_derived_parameters()
        
        self.resolution = resolution #[eV]
        self.dispersion = dispersion #[eV]
        self.collection_angle = collection_angle #[mrad]
        
    def set_derived_parameters(self):
        self.gamma = 1 + self.keV/m0c2 #[au] gamma^2=1-(v/c)^2', label='gamma'
        ################ gamma2=1-v2/c_mat2
        #                      =1-v2*eps_0*eps_r
        #                      =1-eps_r*v2/c_vac2
        ################same as mu in scattering for eps_r=1
        self.beta2 = self.keV * (self.keV + 2 * m0c2)/(self.keV + m0c2)**2 #[au] beta=v/c', label='beta^2' #Note v2/c2=eps_r*v2/c2[vac]
        self.v = c * np.sqrt(self.beta) #[m/s]
        self.T = m0*self.v**2/(2*e) #m0c2/2 * (self.v/c)**2 #[keV]
        
        self.k0 = m0*self.v*self.gamma/(hbar*e) #[1/m]

    def print_parameters(self):
        print('Microscope')
        for each in self.__dict__:
            print(each,': ',self.__dict__[each])
        print()

if __name__ == '__main__':   
    test_scope = Microscope(keV=100)