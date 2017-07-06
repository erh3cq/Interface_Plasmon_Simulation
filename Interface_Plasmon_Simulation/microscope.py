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

class Microscope(HasTraits):
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
    _parent_system = None
    name = Str('User specified', desc='Name of the microscope', label='Name')
    keV = Int(100, desc='Accelerating voltage', label='keV')
    
    beta2 = Float(desc='Relativistic factor, beta=v/c', label='beta^2') #Note v2/c2=eps_r*v2/c2[vac]
    gamma = Float(desc='Relativistic factor in vacuum, gamma^2=1-(v/c)^2', label='gamma')
    ################ gamma2=1-v2/c_mat2
    #                      =1-v2*eps_0*eps_r
    #                      =1-eps_r*v2/c_vac2
    ################same as mu in scattering for eps_r=1
    
    v = Float(desc='Velocity', label='v [m/s]')
    T = Float(desc='Kinetic energy', label='T [keV]')
    k0 = Float(desc='Wave number', label='k0 [1/m]')
    
    resolution = Float(0.05, desc='Energy resolution of the microscope (ie FWHM of the ZLP)', label='Resolution [eV]')
    dispersion = Float(0.05, desc='Energy dispersion of the spectrometer', label='Dispersion [eV]')
    collection_angle = Float(2E-3, desc='Spectrometr collection semi-angle', label='Collection semi-angle [mrad]')
    
    trait_view = View(
            Item('name'),
            'keV',
            Item('gamma', style='readonly', format_str='%.3f'),
            Item('beta', style='readonly', format_str='%.3f'),
            Item('v', style='readonly', format_str='%.3E'),
            Item('T', style='readonly', format_str='%.3E'),
            Item('k0', style='readonly', format_str='%.3E'),
            Group(
                    'resolution',
                    'dispersion',
                    'collection_angle',
                    show_border = True,
                    label='EELS'),
            title = 'Microscope',
            resizable = True)
    
    def __init__(self, resolution=0.05, collection_angle=2E-3, **traits):
        super().__init__(**traits)
#        HasTraits.__init__(self, **traits)
        self.set_derived_parameters()
        
    @on_trait_change('keV')
    def set_derived_parameters(self):
        self.gamma = 1 + self.keV/m0c2 #[au]        
        self.beta = self.keV * (self.keV + 2 * m0c2)/(self.keV + m0c2)**2 #[au]
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
    test_scope.configure_traits()