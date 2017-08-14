# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:38:34 2017

@author: mse334
"""

from traits.api import HasTraits, Float, Str, Instance, on_trait_change, Array
from traitsui.api import Item, View, Label, HGroup
import numpy as np

class vector(HasTraits):
    mag = Float(desc="Lattice parameter")
    hat = Array(np.int, shape=(3,), value=np.array([0,0,0]))
    hat_view = Array(shape=(1,3), value=hat.T)

    def _hat_view_changed(self):
        self.hat = self.hat_view[0].T
        print(self.hat)
        
    traits_view = View(HGroup(Item('mag',show_label=False),Label('[m] * ['),Item('hat_view',show_label=False),Label(']')))
    def __init__(self, **traits):
        super().__init__(**traits)
        if 'hat' in traits: self.hat_view = [traits['hat']]
#        print(self.hat[np.newaxis])

#vector(hat=[0,1,0]).configure_traits()

class EchoBox(HasTraits):
    name = Str('Test')
    a1 = Instance(vector, kw={'hat':[1,0,0]})
    a2 = Instance(vector(hat=[0,1,0]),())
    a3 = Instance(vector, ())
#    a2_mag = Float(desc="Lattice parameter")
#    a2_hat = Array(np.int, shape=(3), value=np.array([0,1,0]))
#    a3_mag = Float(desc="Lattice parameter")
#    a3_hat = Array(np.int, shape=(3), value=np.array([0,0,1]))
    
#    ar = Array(np.int, shape=(3,3), value=np.array([[1,0,0],[0,1,0],[0,0,1]]))
    traits_view = View('name',
                       Item(name='a1', style = 'custom'),
                       Item(name='a2', style = 'custom'),
                       Item(name='a3', style = 'custom'))
    def __init__(self, **traits):
        super().__init__(**traits)
     
box = EchoBox(a3=vector(hat=[0,0,1]))
box.configure_traits()
print(np.zeros((1,3)))