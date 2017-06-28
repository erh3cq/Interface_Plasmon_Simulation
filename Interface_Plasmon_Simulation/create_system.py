# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:04:41 2017

@author: mse334
"""

from traits.api import HasTraits, Float, List, Instance, Str
from traitsui.api import *
from microscope import Microscope
from dimensions import Spectrum_dimensions

#Constants
hbar=6.582E-16#[eV s]

class System(HasTraits):
    name = Str('Main System')
#    _tags = List([microscope, 'dimmension'], label='Tags')
#    _tags2 = List([], label='Tags')
    microscope = Instance(Microscope, desc='Microscope parameters', label='Microscope')
    dimensions = Instance(Spectrum_dimensions, desc='Energy and momentum dimensions', label='Dimmensions')
    
    trait_view = View(
            Item('name', style='readonly'),
            Group(
                    Item(name = 'microscope', style = 'custom', show_label = False, dock = 'tab'),
                    Item(name = 'dimensions', style = 'custom', show_label = False, dock = 'tab'),
                    show_border = True,
                    layout = 'tabbed'),
        title = '',
        resizable = True)
    
    def __init__(self, **traits):
        super().__init__(**traits)
        self.add_microscope(**traits)
        self.dimensions = Spectrum_dimensions(**traits)
        self.dimensions._parent_system = self
        #self.add_trait('_tags2', [self.microscope, 'dimmension', 'test'])

    def add_microscope(self, **traits):
        self.microscope = Microscope(**traits)
        self.microscope._parent_system = self
        
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
    system_test = System(keV=20)
    system_test.configure_traits()