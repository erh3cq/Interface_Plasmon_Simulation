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
                       q_perpendicular=None, E=None,
                       eps=None):
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
    try:
        if eps == 'Drude':
            eps1 = material[0].eps_D #TODO: create heavside function for eps1 and 2
        elif eps == 'Lindhard':
            eps1 = material[0].eps_L
        else:
            eps1 = materials[0].eps
        
        #eps2 = material[1].eps
    except ValueError:
        print('eps has not been set.')
            
    q_parallel = np.sqrt(((E/hbar)/microscope.v)**2+q_perpendicular**2)

    #pre = 2*e**2/(np.pi*ehbar*(microscope.v*q_parallel)**2*eps0) / hbar #[1/(eV m^2)]
    pre = e**2/(np.pi*ehbar*(microscope.v)**2*eps0) / hbar #[1/(eV m^2)]
    
    f_bulk = pre * np.imag(-1/q_parallel**2 * 1/eps1)
    
    return f_bulk

def bulk(microscope, material=None, x=None, x0=0,
                       q_y=None, q_parallel=None, E=None,
                       eps=None):
    """
    First differential scattering cross section d sigma/ dE. 
    """
    
#    x,q_perpendicular,E = np.meshgrid(x,q_perpendicular,E, indexing='ij')
    if q_y is None and q_parallel is None:
        raise Exception('twoSlabParallel requires q_y or q_parallel to be set.')
    elif isinstance(q_y, np.ndarray) is False and isinstance(q_parallel, np.ndarray) is False:
        raise Exception('q_y must be a numpy array.')
    if E is None:
        raise Exception('twoSlabParallel requires E to be set.')
    elif isinstance(E, np.ndarray) is False:
        raise Exception('E must be a numpy array.')
    try:
        if eps == 'Drude':
            eps1 = material[0].eps_D #TODO: create heavside function for eps1 and 2
        elif eps == 'Lindhard':
            eps1 = material[0].eps_L
        elif eps == 'Mermin':
            eps1 = material[0].eps_M
        else:
            eps1 = material[0].eps
        
        #eps2 = material[1].eps
    except ValueError:
        print('eps has not been set.')
            
    gamma2 = 1 - (microscope.v/c)**2 * eps1 #this is acctually the inverse of gamma^2
    if q_parallel is None:
        q_parallel = np.sqrt(((E/hbar)/microscope.v)**2+q_y**2)
    #q2_parallel = ((E/hbar)/microscope.v)**2+q_perpendicular**2
    q_x1 = np.sqrt(q_parallel**2 - ((E/hbar)/c)**2 * eps1)
    

    #pre = 2*e**2/(np.pi*ehbar*(microscope.v*q_parallel)**2*eps0) / hbar #[1/(eV m^2)]
    pre = e**2/(np.pi*ehbar*(microscope.v)**2*eps0) / hbar #[1/(eV m^2)]
    
    f_bulk = pre * np.imag(-1/q_x1 * gamma2/eps1)
    
    return f_bulk

def bulk_q(microscope, material=None, q_perpendicular=None, E=None,
                       eps_tr=None, eps_lon=None):
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
    try:
        if eps_tr == 'Drude':
            eps_tr = material[0].eps_D #TODO: create heavside function for eps1 and 2
        elif eps_tr == 'Lindhard':
            eps_tr = material[0].eps_L
        elif eps_tr == 'Mermin':
            eps_tr = material[0].eps_M
        else:
            eps_tr = material[0].eps
            
        if eps_lon == 'Drude':
            eps_lon = material[0].eps_D #TODO: create heavside function for eps1 and 2
        elif eps_lon == 'Lindhard':
            eps_lon = material[0].eps_L
        elif eps_lon == 'Mermin':
            eps_lon = material[0].eps_M
        else:
            eps_lon = material[0].eps
    except ValueError:
        print('eps has not been set.')
            
    q2 = q_perpendicular**2 + (E/hbar/microscope.v)**2
    pre = e**2/(np.pi*ehbar*(microscope.v)**2*eps0) / hbar #[1/(eV m^2)]
    f_1 = np.imag(1/q2 * (1/eps_tr - 1/eps_lon))
    f_2 = np.imag((microscope.beta2 - 1/eps_tr) * 1/(q2-(E/hbar/c)**2*eps_tr))
    return pre * (f_1+f_2)

def bulk_theta(microscope, material=None, E=None, theta=None):
    """
    """
    E0 = microscope.keV*10E3
    mc2 = m0*c**2/e
    a0 = 4*np.pi*eps0 * hbar**2 / m0 #4*pi*eps0 *hbar^2/(m0*e^2) where hbar is in Js
    
    E_p0 = material.E_p0
    dE = material.dE
    alpha = material.alpha
    
    A = (2 * np.pi**2 * material.na * a0) * 1/E0
    theta_E = E / E0 * (E0 + mc2)/(E0 + 2 * mc2)
    return A * 1/(theta**2 + theta_E**2) * (
            E * dE * E_p0**2 / (
                    (E**2 - E_p0**2 - 4*alpha*E_p0*(theta**2 + theta_E**2))**2 +
                    (E*dE)**2))
    


if __name__ == '__main__':
    from materials import vac, Al
    from Interface_Plasmon_Simulation.microscope import Microscope
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    
    x = np.array([1])#np.arange(-5,0,0.5)*1E-10 #[m]
    E = np.arange(10,20,0.03)
    q_perpendicular = np.arange(100,800,100)*1E6
    
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
    bulk = bulk(microscope, material=materials,
                       x=x, q_perpendicular=q_perpendicular, E=E)
    
    pos_index = -1
    fig, ax = plt.subplots(2, figsize=(20,10), sharex=True)
    fig.suptitle(r'pos: {} [m]'.format(x[pos_index,0,0]))
    img0 = ax[0].imshow(bulk[pos_index,:,:], aspect='auto', origin='lower', norm=LogNorm(), cmap=plt.get_cmap('hot'),
                 extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    ax[0].set_title('Bulk')
    ax[0].set_ylabel(r'$q_{y} [m^-]$')
    #fig.colorbar(mappable=img0, ax=ax[0])
    ax[0].set_aspect('auto')
    
    ax[1].plot(E[0,0,:],bulk[pos_index,:,:].sum(axis=0))
    ax[1].set_xlim(xmin=0)
    plt.autoscale(enable=True, axis='x', tight=True)
    ax[1].set_xlabel('E [eV]')
    #plt.tight_layout(pad=0.001)
    plt.show()