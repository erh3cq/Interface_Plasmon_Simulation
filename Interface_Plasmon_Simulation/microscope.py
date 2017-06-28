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
    _parent_system = None
    name = Str('User specified', desc='Name of the microscope', label='Name')
    keV = Int(100, desc='Accelerating voltage', label='keV')
    gamma = Float(desc='', label='gamma')
    v2_c2 = Float(desc='Speed as a fraction of the speed of light', label='v2/c2')
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
            Item('v2_c2', style='readonly', format_str='%.3f'),
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
        super().__init__()
        self.set_derived_parameters()
        
    @on_trait_change('keV')
    def set_derived_parameters(self):
        self.gamma = 1 + self.keV/m0c2 #[au]        
        self.v2_c2 = self.keV * (self.keV + 2 * m0c2)/(self.keV + m0c2)**2 #[au]
        self.v = c * np.sqrt(self.v2_c2) #[m/s]
        self.T = m0*self.v**2/(2*e) #m0c2/2 * (self.v/c)**2 #[keV]
        
        self.k0 = m0*self.v*self.gamma/(hbar*e) #[1/m]
        
    ##############################
    ##############################
    # Add change catches for keV #
    ##############################
    ##############################

    def print_parameters(self):
        print('Microscope')
        for each in self.__dict__:
            print(each,': ',self.__dict__[each])
        print()

if __name__ == '__main__':   
    test_scope = Microscope(keV=100)
    test_scope.configure_traits()