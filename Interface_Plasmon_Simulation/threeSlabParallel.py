# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 20:34:13 2017

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


def threeSlabParallel_symetric(microscope, material=None, x=None,
                       q_perpendicular=None, E=None, t=1.5E-10):
    """
    Model for an electorn with initial trajectory in the +z direction,
    parallel to a three slab symetric interface that has a normal in the +x direction.
    y is therefore the direction in the plane not parallel to the velocity.
    The electron is in material 1.
    Wang 1996
    
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
            
    q = np.sqrt(((E/hbar)/microscope.v)**2+q_perpendicular**2) #np.sqrt(((E/hbar)/c*microscope.beta)**2+q_perpendicular**2)
    
    #gamma = 1 - eps1 * microscope.beta2
    pre = e**2/(np.pi*ehbar*(microscope.v)**2*eps0) / hbar #[1/(eV m^2)]
    f_n = (eps2+eps1)*(eps2-eps1)* (np.exp(q*t) - np.exp(-1*q*t))
    f_d = (eps2+eps1)**2 * np.exp(q*t) - (eps2-eps1)**2 * np.exp(-1*q*t)
    
    
    f_decay = np.exp(-2*q*np.abs(x-t/2))
    #f_interface_s = pre * np.real(f_decay * ((q_perpendicular/q_x1)**2 * (microscope.v/c)**2 *r_s))
    #f_interface_p = pre * np.real(f_decay * (-1/eps1 * r_p))
    
    return pre*np.imag(1/q * 1/eps1 * f_n/f_d * f_decay)
    

if __name__ == '__main__':
    from Interface_Plasmon_Simulation.materials import vac, Al
    from Interface_Plasmon_Simulation.microscope import Microscope
    import matplotlib.pyplot as plt
    import matplotlib.colors as pltc
    
    t=1.5E-10
    x = np.array([t/2])#np.arange(-5,0,0.5)*1E-10 #[m]
    E = np.arange(0.001,17,0.01)
    q_perpendicular = np.linspace(-1.5,1.5,100,endpoint=True)*1E10#np.arange(0,200,0.05)*1E6
    
    x = x[:,np.newaxis,np.newaxis]
    q_perpendicular = q_perpendicular[np.newaxis,:,np.newaxis]
    E = E[np.newaxis,np.newaxis,:]
    
    materials = [Al, vac]
    for material in materials:
        print(material.name)
        material.set_Ep()
        material.set_eps(E=E)
    
    microscope = Microscope()
    microscope.print_parameters()
    
    interface = threeSlabParallel_symetric(microscope, material=materials,
                       x=x, q_perpendicular=q_perpendicular, E=E, t=t)
    
    norm=pltc.LogNorm()
    fig, ax1 = plt.subplots(1, 1)
    img1 = ax1.imshow(interface[-1,:,:], aspect='auto', norm=norm, origin='lower',
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img1, ax=ax1)
#    plt.colorbar()
    
    plt.show()