# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 18:40:40 2017

@author: erhog
"""

from __future__ import division
import numpy as np

e=1.602E-19 #[C]
c = 3E8 #[m/s]
m0=9.11E-31#[kg]
eps0=8.85E-12#[C^2/N m]
hbar=6.582E-16#[eV s]
ehbar = hbar*e #1.055E-34[J s]

def bulk_NR(microscope, material=None, x=None, x0=0,
                       q_perpendicular=None, E=None):
    """
    First differential scattering cross section d sigma/ dE. 
    """
    
#    x,q_perpendicular,E = np.meshgrid(x,q_perpendicular,E, indexing='ij')
    if q_perpendicular is None:
        raise Exception('twoSlabParallel requires q_perpendicular to be set.')
    elif isinstance(q_perpendicular, np.ndarray) is False:
        raise Exception('q_perpendicular must be a numpy array.')
    if E is None:
        raise Exception('twoSlabParallel requires E to be set.')
    elif isinstance(E, np.ndarray) is False:
        raise Exception('E must be a numpy array.')
    if len(material)!=2:
        raise Exception('material must be a list or array of length 2.')
    else:
        try:
            eps1 = material[0].eps
            eps2 = material[1].eps
        except ValueError:
            print('eps has not been set.')
            
    q_parallel = np.sqrt(((E/hbar)/microscope.v)**2+q_perpendicular**2)

    pre = 2*e**2/(np.pi*ehbar*(microscope.v*q_parallel)**2*eps0) / hbar #[1/(eV m^2)]
    
    f_bulk = pre * np.imag(-1/eps1) #TODO: create heavside function for eps1 and 2
    
    return f_bulk

if __name__ == '__main__':
    from materials import vac, Al
    
    x = np.array([1])#np.arange(-5,0,0.5)*1E-10 #[m]
    E = np.arange(10,50,10)
    q_perpendicular = np.arange(100,800,100)
    
    materials = [Al, vac]
    for material in materials:
        print(material.name)
        material.set_Ep()
        material.set_eps(E=E)
    
    microscope = Microscope()
    microscope.print_parameters()
    bulk_NR(microscope, material=materials,
                       x=x, q_perpendicular=q_perpendicular, E=E)