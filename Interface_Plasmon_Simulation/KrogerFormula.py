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

def kroger(microscope, material=None, x=None, x0=0, t=100E-9,
                       q_perpendicular=None, E=None):
    """
    
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
            
    
    kappa = E/hbar #E/hbar + ky*vy
    q_parallel1 = np.sqrt(q_perpendicular**2 - ((E/hbar)/c)**2 * eps1)#np.sqrt((kappa/microscope.v)**2+q_perpendicular**2)
    q_parallel2 = np.sqrt(q_perpendicular**2 - ((E/hbar)/c)**2 * eps2)
    q_x1 = np.sqrt(q_parallel1**2 + (kappa/microscope.v)**2)
    q_x2 = np.sqrt(q_parallel2**2 + (kappa/microscope.v)**2)
    q_x12 = np.sqrt(q_perpendicular**2 + (kappa/microscope.v)**2 - ((E/hbar)/c)**2 * (eps1 + eps2) )
    
    gamma = 1 - eps2 * microscope.beta2
    
    pre = e**2/(np.pi*ehbar*(microscope.v)**2*eps0) / hbar #[1/(eV m^2)]
    pre_2 = -2*(eps2-eps1)**2/(q_x1**4*q_x2**4)
    
    L_s = q_parallel1*eps2 + q_parallel2*eps1 * np.tanh(q_parallel1*t/2)
    L_as = q_parallel1*eps2 + q_parallel2*eps1 * 1/np.tanh(q_parallel1*t/2)
    B = q_perpendicular*q_x12**2 + eps1*eps2*microscope.beta2*(E/hbar / microscope.v)#*np.sin(theta_n)*np.cos(q_x1)
    
    #f_decay = q_parallel * np.exp(2*q_parallel*np.abs(x-x0))
    #f_interface = pre * np.imag(f_decay * -2/(eps2+eps1))
    #f_begrenzung = pre * np.imag(f_decay * 1/eps1)
    f_bulk = pre * np.imag(1/q_x2**2 * gamma/eps2) * t
    f_interface1 = pre * np.imag(pre_2*(np.sin(kappa/microscope.v*t/2)**2/L_s+np.cos(kappa/microscope.v*t/2)**2/L_as) * B**2/(eps1*eps2))
    
    return f_bulk, f_interface1

    

if __name__ == '__main__':
    from Interface_Plasmon_Simulation.materials import vac, Al
    from Interface_Plasmon_Simulation.microscope import Microscope
    import matplotlib.pyplot as plt
    import matplotlib.colors as pltc
    
    x = np.array([0])*1E-10#np.arange(-5,0,0.5)*1E-10 #[m]
    q_perpendicular = np.arange(0,200,1)*1E6
    E = np.arange(13,17,0.03)
    
    x = x[:,np.newaxis,np.newaxis]
    q_perpendicular = q_perpendicular[np.newaxis,:,np.newaxis]
    E = E[np.newaxis,np.newaxis,:]
    
    materials = [vac, Al]
    for material in materials:
        print(material.name)
        material.set_Ep()
        material.set_eps(E=E)
    
    microscope = Microscope()
    microscope.print_parameters()
    
    
    
    
    bulk, interface1 = kroger(microscope, material=materials,
                       x=x, q_perpendicular=q_perpendicular, E=E)
    print(bulk[0,:,:])


    fig, [ax1, ax2] = plt.subplots(2, 1)
    img1 = ax1.imshow(bulk[0,:,:], aspect='auto', origin='lower',
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img1, ax=ax1)
    
    img2 = ax2.imshow(interface1[0,:,:], aspect='auto', origin='lower',
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img2, ax=ax2)
    
    
    plt.show()