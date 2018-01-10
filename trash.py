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

def kroger(microscope, material=None, x=None, x0=0, t=100E-9, theta_n=0,
                       q_perpendicular=None, E=None):
    """
    
    """
    
    #TODO: microscope.v to vz & vy
    
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
            eps1 = np.conjugate(material[0].eps)
            eps2 = np.conjugate(material[1].eps)
        except ValueError:
            print('eps has not been set.')
            
    
    gamma = 1 - eps2 * microscope.beta2

    ###Angular terms
    theta2 = (q_perpendicular/microscope.k0)**2
    thetaE2 = (E/microscope.keV)#(E/(2*gamma*microscope.T))**2 #(spectrum.E/(hbar*microscope.v*microscope.k0))**2

    lambda2 = theta2 - eps2 * thetaE2 * microscope.beta2
    lambda02 = theta2 - eps1 * thetaE2 * microscope.beta2
    phi2 = lambda2 + thetaE2
    phi02 = lambda02 + thetaE2
    phi2_01 = theta2 + thetaE2 * (1 - (eps2 - eps1) * microscope.beta2)
    
    de = t/2 * E / (hbar * microscope.v)
    tanh = np.tanh(np.sqrt(lambda2/thetaE2) * de)
    L_s = np.sqrt(lambda02) * eps2 + np.sqrt(lambda2) * eps1 * tanh #symetric mode     
    L_as = np.sqrt(lambda02) * eps2 + np.sqrt(lambda2) * eps1 / tanh #anti-symetric mode
    
    
    pre = e**2/(np.pi*ehbar*(microscope.v)**2*eps0) / hbar #[1/(eV m^2)]
    pre_2 = -2 * theta2 * (eps2 - eps1)**2 / (microscope.k0 * (phi02 * phi2)**4)#2)#
        
    f_bulk = pre * np.imag(t * gamma / (eps2 * phi2))
    
    f_surface = pre * np.imag(pre_2 * phi2_01**2 / (eps2*eps1)
                               * (np.sin(de)**2/L_s + np.cos(de)**2/L_as))
    f_guidedLight1 = pre * np.imag(pre_2 * microscope.beta2 * np.sqrt(lambda02*thetaE2) * phi2_01 / eps1 *
                                (1/L_s - 1/L_as) * np.sin(2 * de))
        
    f_guidedLight2 = pre * np.imag(pre_2 * -1* microscope.beta2**2 * np.sqrt(lambda02*lambda2) * thetaE2 *
                                (np.cos(de)**2 * tanh/L_s + np.sin(de)**2/(L_as * tanh)))
    
    
    
    
    
    
    
    
    
    return f_bulk, f_surface, f_guidedLight1, f_guidedLight2#f_bulk, f_interface1, f_interface2

    

if __name__ == '__main__':
    from Interface_Plasmon_Simulation.materials import vac, Al
    from Interface_Plasmon_Simulation.microscope import Microscope
    import matplotlib.pyplot as plt
    import matplotlib.colors as pltc
    
    x = np.array([0])*1E-10#np.arange(-5,0,0.5)*1E-10 #[m]
    q_perpendicular = np.arange(0,200,1)*1E6#np.linspace(-1.3,1.3,20)*1E10#
    E = np.arange(1E-3,25,0.03)
    
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
    
    
    norm = None
    bulk, surface, guidedLight1, guidedLight2 = kroger(microscope, material=materials, t=30E-9, theta_n=0,
                       x=x, q_perpendicular=q_perpendicular, E=E)


    fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(5, 1)
    img1 = ax1.imshow(bulk[0,:,:], aspect='auto', origin='lower', norm=pltc.LogNorm(),
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img1, ax=ax1)
    
    img2 = ax2.imshow(surface[0,:,:], aspect='auto', origin='lower', norm=pltc.LogNorm(),
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img2, ax=ax2)
    
    img3 = ax3.imshow(guidedLight1[0,:,:], aspect='auto', origin='lower', norm=pltc.LogNorm(),
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img3, ax=ax3)
    
    img4 = ax4.imshow(guidedLight2[0,:,:], aspect='auto', origin='lower', norm=norm,
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img4, ax=ax4)
    
    img5 = ax5.imshow(bulk[0,:,:]+surface[0,:,:]+guidedLight1[0,:,:]+guidedLight2[0,:,:], aspect='auto', origin='lower', norm=pltc.LogNorm(),
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img5, ax=ax5)
    
    fig.suptitle('Egerton')
    plt.show()