# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:38:34 2017

@author: mse334
"""

from traits.api import HasTraits, HasPrivateTraits, Float, Enum, Str, Instance, on_trait_change, Any, Generic,Undefined, Trait
from traitsui.api import *

class Microscope(HasTraits):
#    parent_system = Instance(Undefined)#Generic
    _parent_system = None
    inp = Float(0)
    out = Float(0)
        
#    @on_trait_change('parent_system')
#    def _parent_sys_changed(self):
#        self.__parent_sys = self
#        print('parent system set to: ', self.parent_system.name)

    @on_trait_change('inp')
    def _parent_sysN_changed(self):
        print('Input(scope) changed to: ', self.inp)
        self.out = self.inp
        
#   WORKS
#    def _parent_sysN_changed(self):
#        print('parent system set to: ', self.parent_system.name)
#        self.out = self.inp    

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.add_trait('name',Enum('Titan','Themis','Nion'))

    trait_view = View(Item(name = 'name'),
#                      Item(name = 'parent_system'),
                      Group(Item(name = 'inp'),
                            Item(name = 'out'),
                            label = '%s parameters'%inp,
                            show_border = True))
        

class Dimensions(HasTraits):
    def __init__(self):
        self.number = 'two'

class System(HasTraits):
    name = Str('Main System')
    scope = Instance(Microscope, (), desc='Microscope parameters', label='Microscope')
    dimmension = Instance(Dimensions, desc='Energy and momentum dimensions', label='Dimmensions')
    
    inp = Float(0, label='Input')
    output = Float(0)
    @on_trait_change('inp')
    def _inp_changed(self):
        self.output = self.inp
        print('Input (system) changed to',self.inp)
    
    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        self.initialized = 'Yes'
        self.add_microscope()
        self.axis = Dimensions()
        
    trait_view = View(
            Item('name', style='readonly'),
            Group(
                    Item(name = 'scope', style = 'custom', show_label = False, dock = 'tab'),
                    Group(
                            Item(name = 'inp'),
                            Item(name = 'output'),
                            label = 'System parameters',
                            dock = 'tab'),
                    show_border = True,
                    layout = 'tabbed'),
        title = '',
        resizable = True)
        
        
    def add_microscope(self):
#        self.scope = Microscope()
#        self.scope.parent_system.klass = self
#        self.scope.parent_system = self
        Microscope.add_class_trait('parent_system',Instance(self))
        self.scope.parent_system = self
        self.scope._parent_system = self
        print('Added scope with value of:',self.scope.parent_system)
        
           
tester = System()
#print(tester.trait('parent_system').default_kind) 
tester.configure_traits()

print('\n')

class EchoBox(HasPrivateTraits):
    pre_inp = Float(label='Pre-input')
    _inp1 = Float(1)
    _inp2 = Float(2)
    output = Float()
    trait_view = View(Group(Item(name = '_inp1'),
                            Item(name = '_inp2'),
                             Item(name = 'output'),
                             label = 'Section',
                             show_border = True))
#    def _inp_changed(self):
#        self.output = self._inp
#        print('changed')
    @on_trait_change(['_inp1','_inp2'])
    def _inp_changed(self):
        self.output = self._inp1+self._inp2

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

box.configure_traits()
    