# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 12:43:02 2017

@author: mse334
"""

import numpy as np
from traits.api import HasTraits, Float, Int, Generic, Instance, on_trait_change
from traitsui.api import *

#Constants
hbar=6.582E-16#[eV s]


class Spectrum_dimensions(HasTraits):
    """
    Class containint the dimensions of the spectrum.
    
    Adjustable parameters
    ----------------------
    E_min : Float
        Minimum energy loss of spectrum.
        Units in eV.
    E_max : Float
        Maximum energy loss of spectrum.
        Units in eV.
    E_N : Int
        Number of points to calcualte in the energy range.
        Units in eV.
    q_perpendicular_max : Float
        Maximum momentum to calculate up to, that is transfered perpendicular to the incident beam.
        Units in 1/m.
    q_perpendicular_N : Int
        Number of points to calcualte in the q_perpendicular range.
        Units in 1/m.
        
    Derived parameters
    -------------------
    E : Array
        Array from E_min to E_max consisting of E_N points.
        Units in eV.
    q_perpendicular : Array
        Array from \+- q_perpendicular_max with q_perpendicular_N points.
        Units in 1/m.
    """
    _parent_system = Generic
    E_min = Float(5, desc='Minimum energy loss of spectrum', label='Minimum')
    E_max = Float(20, desc='Maximum energy loss of spectrum', label='Maximum')
    E_N = Int(400, desc='Number of points to calcualte in the energy range', label='N')
    q_perpendicular_max = Float(1E10, desc='Maximum momentum to calculate up to, that is transfered perpendicular to the incident beam.', label='Maximum')
    q_perpendicular_N = Int(800, desc='Number of points to calcualte in the q_perpendicular range', label='N')

    trait_view = View(
            '',
            Group(
                    Item('E_min'),
                    Item('E_max'),
                    Item('E_N'),
                    label = 'Energy',
                    show_border = True),
            Group(
                    Item('q_perpendicular_max'),
                    Item('q_perpendicular_N'),
                    label='Momentum',
                    show_border = True),
            dock = 'vertical',
            title = 'Spectrum dimensions',
            resizable = True)

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.E = np.linspace(self.E_min, self.E_max, self.E_N)
        self.q_perpendicular = np.linspace(-1.0*self.q_perpendicular_max, self.q_perpendicular_max, self.q_perpendicular_N)[:,None]#microscope.k0*microscope.collection_angle
        self.q_parallel = None
    
    @on_trait_change(['E_min','E_max', 'E_N'])
    def _energy_range_changed(self):
        self.E = np.linspace(self.E_min, self.E_max, self.E_N)
#        print(self.E)
        
    @on_trait_change(['q_perpendicular_max', 'q_perpendicular_N'])
    def _q_perpendicular_range_changed(self):
        self.q_perpendicular = np.linspace(-1.0*self.q_perpendicular_max, self.q_perpendicular_max, self.q_perpendicular_N)[:,None]
#        print('q_perpendicular updated')

    @on_trait_change('_parent_system')
    def _parent_system_changed(self):
        print('Parent change',self._parent_system)
        
#    @property
#    def q_parallel(self):
##        print('getting q_parallel')
#        return self.__q_parallel
#    @q_parallel.setter
#    def q_parallel(self,value):
#        if value is None:
#            if self.__parent.microscope:
##                print(self.__parent.microscope.v)
##                print('Set q_parallel to ',self.E/(hbar * self.__parent.microscope.v),' based on microscope parameters')
#                self.__q_parallel = self.E/(hbar * self.__parent.microscope.v)
#            else:
##                print('Microscope parameters have not been set.  Therefore, q_parallel = None')
#                self.__q_parallel = None
#        else:
##            print('Set q_parallel to specified value of ',value)
#            self.__q_parallel = value

if __name__ == '__main__':
    dimension_test = Spectrum_dimensions()
    dimension_test.configure_traits()