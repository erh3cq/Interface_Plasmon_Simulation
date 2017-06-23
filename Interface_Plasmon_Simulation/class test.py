# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:38:34 2017

@author: mse334
"""

from traits.api import HasTraits, Float, List, Str, Instance, on_trait_change, Property
from traitsui.api import *

class Microscope(HasTraits):
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
        

class Dimensions(HasTraits):
    def __init__(self):
        self.number = 'two'

class System(HasTraits):
    name = Str('Main system')
    microscope = Instance(Microscope, desc='Microscope parameters', label='Microscope')
    dimmension = Instance(Dimensions, desc='Energy and momentum dimensions', label='Dimmensions')
    def __init__(self):
        self.initialized = 'Yes'
        self.add_microscope()
        self.axis = Dimensions()
        
    def add_microscope(self):
        self.microscope = Microscope()
#        self.microscope.parent_sys = self
        print(self.microscope.parent_sys)
        Microscope.add_class_trait('parent_sys',self)
        print(self.microscope.parent_sys)

            
tester = System()
#tester.configure_traits()

print('\n')

class EchoBox(HasTraits):
    pre_inp = Float
    inp = None
    output = Float(0)
    def _inp_changed(self):
        self.output = self.inp
        print('changed')
    def add_inp(self):
        EchoBox.add_class_trait('inp',Float(1., depends_on='pre_inp'))
        
box = EchoBox()
print(box.inp)
box.add_inp()
print(box.inp)
box.inp = 2
print(box.inp)
#box.configure_traits()
    