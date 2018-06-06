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


def twoSlabParallel(microscope, material=None, x=None, x0=0,
                       q_y=None, q_parallel=None, E=None, NR=False):
    """
    Model for an electorn with initial trajectory in the +z direction,
    parallel to an interface that has a normal in the +x direction.
    y is therefore the direction in the plane not parallel to the velocity.
    The electron is in material 1.
    Wang 1996
    
    q_xi: Momentum transfer of a SP in the direction normal to an interface.
    q_parallel: Momentum transfer of a SP parallel to the electrons trajectory.
    q_y: Momentum transfer perpendicular
                    to the interface normal and the beam (ie. q_x x q_parallel)
    """
    
#    x,q_perpendicular,E = np.meshgrid(x,q_perpendicular,E, indexing='ij')
    if NR is True:
        print('NR')
        microscope.beta2 = 0
    if q_y is None and q_parallel is None:
        raise Exception('twoSlabParallel requires q_y or q_parallel to be set.')
    elif isinstance(q_y, np.ndarray) is False and isinstance(q_parallel, np.ndarray) is False:
        raise Exception('q_y must be a numpy array.')
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
    
    gamma2 = 1 - eps1 * microscope.beta2
    
    if q_parallel is None:
        q_parallel = np.sqrt(((E/hbar)/microscope.v)**2+q_y**2)
    q_x1 = np.sqrt(q_parallel**2 - ((E/hbar)/microscope.v)**2 * eps1 * microscope.beta2)
    q_x2 = np.sqrt(q_parallel**2 - ((E/hbar)/microscope.v)**2 * eps2 * microscope.beta2)
    pre = e**2/(2*np.pi**2*ehbar*(microscope.v)**2*eps0) / hbar #[1/(eV m^2)]
    
    #f_1 = 2*(eps2 - eps1)*q_x1 / (eps1*(q_x1+q_x2)*(eps2*q_x1 + eps1*q_x2))
    #f_2 = -1*gamma/(eps1*q_x1) * (q_x1 - q_x2)/(q_x1 + q_x2)
    
    f_decay = np.real(np.exp(-2*q_x1*np.abs(x-x0)))
    f_pre = 1/(q_x1 + q_x2)
    f_1 = pre * np.imag(-2 * (q_x1 + q_x2)/(q_x1*eps2 + q_x2*eps1) * f_pre)* f_decay

    f_2 = pre * np.imag(2/eps1 * f_pre)* f_decay
    f_3 = pre * np.imag(-1*gamma2/eps1 * (q_x1 - q_x2)/q_x1 * f_pre)* f_decay
    
    return f_1, f_2, f_3

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
    
    interface, begrenzung = twoSlabParallel_NR(microscope, material=materials,
                       x=x, q_perpendicular=q_perpendicular, E=E)
    
    fig, [ax1, ax2] = plt.subplots(2, 1)
    img1 = ax1.imshow(interface[-1,:,:], aspect='auto', norm=pltc.LogNorm(), origin='lower',
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img1, ax=ax1)
#    plt.colorbar()
    
    img2 = ax2.imshow(begrenzung[-1,:,:], aspect='auto', origin='lower',
                 extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img2, ax=ax2)
    
    interface_s, interface_p = twoSlabParallel(microscope, material=materials,
                       x=x, q_perpendicular=q_perpendicular, E=E)
    
    fig, [ax1, ax2] = plt.subplots(2, 1)
    img1 = ax1.imshow(interface_s[-1,:,:], aspect='auto', origin='lower',
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img1, ax=ax1)
#    plt.colorbar()
    
    img2 = ax2.imshow(interface_p[-1,:,:], aspect='auto', origin='lower',
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