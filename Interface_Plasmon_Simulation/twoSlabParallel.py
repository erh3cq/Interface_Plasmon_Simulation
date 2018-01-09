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

def twoSlabParallel_NR(microscope, material=None, x=None, x0=0,
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
            
    q_parallel = np.sqrt(((E/hbar)/microscope.v)**2+q_perpendicular**2)

    pre = 2*e**2/(np.pi*ehbar*(microscope.v*q_parallel)**2*eps0) / hbar #[1/(eV m^2)]
    
    
    f_decay = q_parallel * np.exp(2*q_parallel*np.abs(x-x0))
    f_interface = pre * np.imag(f_decay * -2/(eps2+eps1))
    f_begrenzung = pre * np.imag(f_decay * 1/eps1)
    
    return f_interface, f_begrenzung

def twoSlabParallel(microscope, material=None, x=None, x0=0,
                       q_perpendicular=None, E=None):
    """
    Model for an electorn with initial trajectory in the +z direction,
    parallel to an interface that has a normal in the +x direction.
    y is therefore the direction in the plane not parallel to the velocity.
    The electron is in material 1.
    Garcia de 
    
    q_xi: Momentum transfer in the direction parallel to the beam in material i.
    q_parallel: Momentum transfer parallel to the interface normal.
    q_perpendicular: Momentum transfer perpendicular
                    to the interface normal and the beam (ie. q_x x q_parallel)
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
            
    q_parallel = np.sqrt(((E/hbar)/microscope.v)**2+q_perpendicular**2) #np.sqrt(((E/hbar)/c*microscope.beta)**2+q_perpendicular**2)
    q_x1 = np.sqrt(((E/hbar)/c)**2 * eps1 - q_parallel**2)
    q_x2 = np.sqrt(((E/hbar)/c)**2 * eps2 - q_parallel**2)
    pre = 2*e**2/(np.pi*ehbar*(microscope.v*q_parallel)**2*eps0) / hbar #[1/(eV m^2)]
    r_p = (eps2*q_x1 - eps1*q_x2)/(eps2*q_x1 + eps1*q_x2)
    r_s = (q_x1 - q_x2)/(q_x1 + q_x2)
    
    
    f_decay = q_x1 * np.exp(2j*q_x1*np.abs(x-x0))
    f_interface_s = pre * np.real(f_decay * ((q_perpendicular/q_x1)**2 * (microscope.v/c)**2 *r_s))
    f_interface_p = pre * np.real(f_decay * (-1/eps1 * r_p))
    
    return f_interface_s, f_interface_p

def twoSlabParallel_2(microscope, material=None, x=None, x0=0,
                       q_perpendicular=None, E=None):
    """
    Model for an electorn with initial trajectory in the +z direction,
    parallel to an interface that has a normal in the +x direction.
    y is therefore the direction in the plane not parallel to the velocity.
    The electron is in material 1.
    Wang 1996
    
    q_xi: Momentum transfer in the direction parallel to the beam in material i.
    q_parallel: Momentum transfer parallel to the interface normal.
    q_perpendicular: Momentum transfer perpendicular
                    to the interface normal and the beam (ie. q_x x q_parallel)
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
            
    q_parallel = np.sqrt(((E/hbar)/microscope.v)**2+q_perpendicular**2) #np.sqrt(((E/hbar)/c*microscope.beta)**2+q_perpendicular**2)
    q_x1 = np.sqrt(q_parallel**2 - ((E/hbar)/c)**2 * eps1)
    q_x2 = np.sqrt(q_parallel**2 - ((E/hbar)/c)**2 * eps2)
    gamma = 1 - eps1 * microscope.beta2
    pre = e**2/(np.pi*ehbar*(microscope.v)**2*eps0) / hbar #[1/(eV m^2)]
    f_1 = 2*(eps2 - eps1)*q_x1 / (eps1*(q_x1+q_x2)*(eps2*q_x1 + eps1*q_x2))
    f_2 = -1*gamma/(eps1*q_x1) * (q_x1 - q_x2)/(q_x1 + q_x2)
    
    
    f_decay = np.exp(-2*q_x1*np.abs(x-x0))
    #f_interface_s = pre * np.real(f_decay * ((q_perpendicular/q_x1)**2 * (microscope.v/c)**2 *r_s))
    #f_interface_p = pre * np.real(f_decay * (-1/eps1 * r_p))
    
    return pre*np.imag(f_1*f_decay), pre*np.imag(f_2*f_decay)

def dispersion(material=None, E=None):
    eps1 = material[0].eps
    eps2 = material[1].eps
    disp = E/(hbar*c)*np.sqrt(eps1*eps2/(eps1+eps2))
    return np.real(disp)
    

if __name__ == '__main__':
    from Interface_Plasmon_Simulation.materials import vac, Al
    from Interface_Plasmon_Simulation.microscope import Microscope
    import matplotlib.pyplot as plt
    import matplotlib.colors as pltc
    
    
    x = np.array([0])#np.arange(-5,0,0.5)*1E-10 #[m]
    E = np.arange(0.001,17,0.01)
    q_perpendicular = np.arange(0,200,0.05)*1E6
    
    materials = [Al, vac]
    for material in materials:
        print(material.name)
        material.set_Ep()
        material.set_eps(E=E)
    
    microscope = Microscope()
    microscope.print_parameters()
    
    interface, begrenzung = twoSlabParallel_NR(microscope, material=materials,
                       x=x, q_perpendicular=q_perpendicular, E=E)
    
    fig, [ax1, ax2] = plt.subplots(2, 1)
    img1 = ax1.imshow(interface[:,-1,:], aspect='auto', norm=pltc.LogNorm(), origin='lower',
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img1, ax=ax1)
#    plt.colorbar()
    
    img2 = ax2.imshow(begrenzung[:,-1,:], aspect='auto', origin='lower',
                 extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img2, ax=ax2)
    
    interface_s, interface_p = twoSlabParallel(microscope, material=materials,
                       x=x, q_perpendicular=q_perpendicular, E=E)
    
    fig, [ax1, ax2] = plt.subplots(2, 1)
    img1 = ax1.imshow(interface_s[:,-1,:], aspect='auto', origin='lower',
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img1, ax=ax1)
#    plt.colorbar()
    
    img2 = ax2.imshow(interface_p[:,-1,:], aspect='auto', origin='lower',
                 extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img2, ax=ax2)
    ax2.plot(q_perpendicular*c*hbar, q_perpendicular, color='red')
    ax2.set_xlim(0,np.amax(E))
    
    fig2 = plt.figure('dispersion')
    plt.plot(E,dispersion(material=materials, E=E))
    plt.plot(q_perpendicular*c*hbar, q_perpendicular, color='red')
    plt.xlim(0,np.amax(E))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{||} [m^-]$')
    
    plt.show()