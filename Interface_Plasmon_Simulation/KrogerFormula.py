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
            
    vx = microscope.v * np.cos(np.pi/180*theta_n)
    vy = microscope.v * np.sin(np.pi/180*theta_n)
    kappa = E/hbar #E/hbar + ky*vy
    print(kappa)
    print(vx)
    print(kappa/vx)
    cosq_x2 = np.cos(0)
    
    gamma = 1 - eps2 * microscope.beta2
    
    q_parallel1 = np.sqrt(q_perpendicular**2 - ((E/hbar)/c)**2 * eps1)
    q_parallel2 = np.sqrt(q_perpendicular**2 - ((E/hbar)/c)**2 * eps2)
    q_x1 = np.sqrt(q_parallel1**2 + (kappa/vx)**2)
    q_x2 = np.sqrt(q_parallel2**2 + (kappa/vx)**2)
    q_x12 = np.sqrt(q_perpendicular**2 + (kappa/vx)**2 - ((E/hbar)/c)**2 * (eps1 + eps2) )
    
    pre = e**2/(np.pi*ehbar*(microscope.v)**2*eps0) / hbar #[1/(eV m^2)]
    pre_2 = -2*(eps2-eps1)**2/(q_x1**4*q_x2**4)
    
    L_s = q_parallel1*eps2 + q_parallel2*eps1 * np.tanh(q_parallel1*t/2)
    L_as = q_parallel1*eps2 + q_parallel2*eps1 * 1/np.tanh(q_parallel1*t/2)
    A = (q_perpendicular - kappa/vx * np.tan(np.pi/180*theta_n) * cosq_x2) * \
        E/hbar/microscope.v*microscope.beta2*np.cos(np.pi/180*theta_n)
    B = q_perpendicular*q_x12**2 + eps1*eps2*microscope.beta2**2 * (E/hbar / microscope.v)**3 * np.sin(np.pi/180*theta_n)*cosq_x2
    
    
    f_bulk = pre * np.imag(1/q_x2**2 * gamma/eps2 * t)
    f_interface = pre * np.imag(pre_2 * (
            np.sin(kappa/vx * t/2)**2/L_s + np.cos(kappa/vx * t/2)**2/L_as) * \
            B**2/(eps1*eps2))
    f_guidedLight1 = pre * np.imag(-1*pre_2 * (
            np.cos(kappa/vx * t/2)**2/L_s * np.tanh(q_parallel2 * t/2) +
                                        np.sin(kappa/vx * t/2)**2/L_as * 1/np.tanh(q_parallel2  *t/2)) *
            A**2*q_parallel1*q_parallel2)
    f_guidedLight2 = pre * np.imag(pre_2 * (
            (1/L_s - 1/L_as) * q_parallel1/eps1 * A*B * np.sin(t * kappa/vx)))
    
    
    return f_bulk, f_interface, f_guidedLight1, f_guidedLight2

    

if __name__ == '__main__':
    from Interface_Plasmon_Simulation.materials import vac, Al
    from Interface_Plasmon_Simulation.microscope import Microscope
    import matplotlib.pyplot as plt
    import matplotlib.colors as pltc
    
    x = np.array([0])*1E-10#np.arange(-5,0,0.5)*1E-10 #[m]
    q_perpendicular = np.linspace(0,200,20)*1E6#np.arange(0,200,1)*1E6#np.linspace(-1.3,1.3,20)*1E10#
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
    
    
    t=10E-9
    theta_n=0
    norm = pltc.LogNorm()
    bulk, interface, guidedLight1, guidedLight2 = kroger(microscope, material=materials, t=t, theta_n=theta_n,
                       x=x, q_perpendicular=q_perpendicular, E=E)


    fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(5, 1)
    fig.suptitle(r't: {}[nm],  $\theta_n$: {}$^o$'.format(t/1E-9,theta_n))
    
    img1 = ax1.imshow(bulk[0,:,:], aspect='auto', origin='lower', norm=norm,
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img1, ax=ax1)
    
    img2 = ax2.imshow(interface[0,:,:], aspect='auto', origin='lower', norm=norm,
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img2, ax=ax2)
    
    img3 = ax3.imshow(guidedLight1[0,:,:], aspect='auto', origin='lower', norm=norm,
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img3, ax=ax3)
    
    img4 = ax4.imshow(guidedLight2[0,:,:], aspect='auto', origin='lower', norm=norm,
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img4, ax=ax4)
    
    img5 = ax5.imshow(bulk[0,:,:]+interface[0,:,:]+guidedLight1[0,:,:]+guidedLight2[0,:,:],
                      aspect='auto', origin='lower', norm=norm,
                     extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
    plt.xlabel('E [eV]')
    plt.ylabel(r'$q_{y} [m^-]$')
    fig.colorbar(img5, ax=ax5)
    
    plt.show()