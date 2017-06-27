# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 13:31:02 2017

@author: mse334
"""
import numpy as np
from traits.api import HasTraits, Float, List, Instance, on_trait_change
from traitsui.api import *
        
class create_sample():
    def __init__(self, slabs, spectrum_dimmensions, interface_angle=0):
        self.slabs = slabs
        self.interface_angle = interface_angle
        self.set_slab_positions()
        self.set_Ep_and_eps(spectrum_dimmensions.E ,spectrum_dimmensions.q_perpendicular)
    def set_slab_positions(self):
        for index, slab in enumerate(self.slabs):
            if index==0:
                slab.position = -np.inf
            elif index==1:
                slab.position = 0
            else:
                slab.position = slabs[index-1].position + slabs[index-1].width
    def set_Ep_and_eps(self, E ,q_perpendicular):
        for slab in self.slabs:
            slab.material.set_Ep(q=q_perpendicular)
            slab.material.set_eps(E=E)
    def print_details(self, slab='All'):
        def _info(index):
            print('_________________')
            print('Slab ',index,':')
            print('position: %.2f [nm]   , width: %.2f [nm]'%(self.slabs[index].position/10**-9, self.slabs[index].width/10**-9))
            print('Ep_0:   %.2f [eV]   , dE:    %.2f [eV]'%(self.slabs[index].material.Ep_0, self.slabs[index].material.dE))
            print('q_min:    %.2f [1/nm] , q_c:   %.2f [1/nm]'%(self.slabs[index].q_min/10**9, self.slabs[index].q_c/10**9))    
        if slab == 'All':
            for index, slab in enumerate(self.slabs):
                _info(index)
        else:
            _info(slab)
            
            
#class slab_system():
#    def __init__(self, slabs, z, energy ,q_perpendicular, interface_angle=2E-3):
#        self.slabs = slabs
#        self.interface_angle = interface_angle
#        self.set_slab_positions()
#        self.set_Ep_and_eps(energy ,q_perpendicular)
#    def set_slab_positions(self):
#        for index, slab in enumerate(self.slabs):
#            if index==0:
#                slab.position = -np.inf - z * self.interface_angle
#            elif index==1:
#                slab.position = 0 - z * self.interface_angle
#            else:
#                slab.position = slabs[index-1].position + slabs[index-1].width# - z * self.interface_angle
#    def set_Ep_and_eps(self, energy ,q_perpendicular):
#        for slab in self.slabs:
#            slab.material.set_Ep(q=q_perpendicular)
#            slab.material.set_eps(E=energy)
#    def print_details(self, slab='All'):
#        def _info(index):
#            print('_________________')
#            print('Slab ',index,':')
#            print('position: %.2f [nm]   , width: %.2f [nm]'%(self.slabs[index].position[0]/10**-9, self.slabs[index].width/10**-9))
#            print('Ep_0:   %.2f [eV]   , dE:    %.2f [eV]'%(self.slabs[index].material.Ep_0, self.slabs[index].material.dE))
#            print('q_min:    %.2f [1/nm] , q_c:   %.2f [1/nm]'%(self.slabs[index].q_min/10**9, self.slabs[index].q_c/10**9))    
#        if slab == 'All':
#            for index, slab in enumerate(self.slabs):
#                _info(index)
#        else:
#            _info(slab)