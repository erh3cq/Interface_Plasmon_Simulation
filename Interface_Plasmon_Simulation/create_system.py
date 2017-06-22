# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:04:41 2017

@author: mse334
"""

from traits.api import HasTraits, Float, List, Instance
from traitsui.api import *
from microscope import microscope
from create_slabs import spectrum_dimmensions

#Constants
hbar=6.582E-16#[eV s]

class system(HasTraits):
    
#    _tags2 = List([], label='Tags')
    test = Float(2.0, label='Test')
    microscope = Instance(microscope, desc='Microscope parameters', label='Microscope')
    dimmension = Instance(spectrum_dimmensions, desc='Energy and momentum dimensions', label='Dimmensions')
    #_tags = List([microscope, 'dimmension'], label='Tags')
    def __init__(self):
        self.add_microscope()
        self.dimmension = spectrum_dimmensions(self)
        spectrum_dimmensions.add_class_trait('parent_system',self)
        #self.add_trait('_tags2', [self.microscope, 'dimmension', 'test'])
        #self.dimmension = Instance(spectrum_dimmensions(self), desc='Energy and momentum dimensions', label='Dimmensions')

    def add_microscope(self):
        self.microscope = microscope()
        microscope.parent_system = self
#        microscope.add_class_trait('parent_system',self)
        
    def add_tag(self):
        return
    def print_tags(self):
        print('System')
#        tag_values = [each for each in self.__dict__ if each in self.tags]
        for tag in self.__dict__:
            print('|-',tag)
            for tag2 in tag.__dict__:
                print('  |-',tag2)


if __name__ == '__main__':
    System = system()
    print('Print E_min',System.dimmension.E_min)
    System.dimmension.E_min = 10.0
    #print(system._tags2[2])
    
    System.configure_traits()