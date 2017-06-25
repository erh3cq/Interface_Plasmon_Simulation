# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:38:34 2017

@author: mse334
"""

from traits.api import HasTraits, HasPrivateTraits, Float, List, Str, Instance, on_trait_change, Any, Generic
from traitsui.api import *

class Microscope(HasTraits):
    parent_sys = Generic
    parent_sysN = Float(0)
    out = Float(0)
        
    @on_trait_change('parent_sys')
    def _parent_sys_changed(self):
        self.__parent_sys = self
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
#        Microscope.add_class_trait('parent_sys',self)
        print(self.microscope.parent_sys)

            
tester = System()
#tester.configure_traits()

print('\n')

class EchoBox(HasPrivateTraits):
    pre_inp = Float(label='Pre-input')
    _inp = Generic
    output = Float(0)
    trait_view = View(Group(Item(name = 'pre_inp'),
                             Item(name = 'output'),
                             label = 'Section',
                             show_border = True))
#    def _inp_changed(self):
#        self.output = self._inp
#        print('changed')
    @on_trait_change(['_inp','pre_inp'])
    def _inp_changed(self):
        self.output = self.pre_inp
        print('changed to',self.pre_inp)

    def add_inp(self):
        if self._inp is not self.trait('_inp').is_trait_type( Float ):
            print('not Float conditional')
            print('1:',self.trait('_inp').trait_type)
#            self._inp = Float(0.0)
            self._inp = 0.0
            print('2:',self.trait('_inp').trait_type)
        else:
            print('Float Conditional')
            self._inp = self._inp+1
#        EchoBox.add_class_trait('_inp',Float(1., depends_on='pre_inp'))
        
box = EchoBox()
print(box._inp)
print('Is float: ',box.trait('_inp').is_trait_type( Float ) )
print('Is any: ',box.trait('_inp').is_trait_type( Any ) )
print()

box.add_inp()
print('3',box._inp)
box.trait('_inp').Trait_type=Float
print(box.trait('_inp').Trait_type )
print('4',box._inp)
print('Is float: ',box.trait('_inp').is_trait_type( Float ) )
print('Is any: ',box.trait('_inp').is_trait_type( Any ) )
print()


box.add_inp()
#box._inp = 2
print('3',box._inp)
print('Is float: ',box.trait('_inp').is_trait_type( Float ) )
print('Is any: ',box.trait('_inp').is_trait_type( Any ) )
box.configure_traits()
    