# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:38:34 2017

@author: mse334
"""

from traits.api import HasTraits, Float, List, Str, Instance, on_trait_change
from traitsui.api import *

class microscope(HasTraits):
    parent_sys = None
    parent_sysN = Float(0)
    out = Float(0)
#    @property
#    def parent_sys(self):
#        return self.parent_sys
#    @parent_sys.setter
#    def parent_sys(self, parent):
#        self.__parent_sys = parent
#        print('parent system set to: ',parent)
        
    @on_trait_change('parent_sysN')
    def _parent_sysN_changed(self):
#        self.__parent_sys = parent
        print('parent system set to: ', self.parent_sys.name)
        self.out = self.parent_sysN
        
#   WORKS
#    def _parent_sysN_changed(self):
#        print('parent system set to: ', self.parent_sys.name)
#        self.out = self.parent_sysN    

    def __init__(self):
        self.name = Str('Titan')
        

class dimensions(HasTraits):
    def __init__(self):
        self.number = 'two'

class system(HasTraits):
    name = Str('Main system')
    microscope = Instance(microscope, desc='Microscope parameters', label='Microscope')
    dimmension = Instance(dimensions, desc='Energy and momentum dimensions', label='Dimmensions')
    def __init__(self):
        self.initialized = 'Yes'
        self.add_microscope()
        self.axis = dimensions()
        
    def add_microscope(self):
        self.microscope = microscope()
#        self.microscope.parent_sys = self
        microscope.add_class_trait('parent_sys',self)

            
tester = system()
tester.configure_traits()