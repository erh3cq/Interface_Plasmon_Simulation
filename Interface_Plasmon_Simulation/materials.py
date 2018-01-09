# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 17:11:43 2016

@author: Eric Hoglund
"""

from __future__ import division
import numpy as np
import scipy as sp

hbar=6.582E-16#[eV s]   1.055E-34[J s]
e=1.602E-19 #[C]
m0=9.11E-31#[kg]
eps0=8.85E-12#[C^2/N m]


class material():
    #TODO: add effective mass
    def __init__(self, name='General Material', a=1E-9, dz=1,
                 lattice = [[1, 0, 0],[0, 1, 0],[0, 0, 1]],
                 basis = [[0,0,0]],
                 valance=3, dE=10, me=1):
        self.name=name
        if isinstance(a, list):
            self.a = np.array(a)
        else:
            self.a = a*np.ones(3) #[m]
        
        self.dz = dz #density adjustment factor (ex. on average 1 of 4 atoms in FCC is missing then dz=0.75)
        self.valance = valance #valance per basis
        self.dE = dE
#        self.damp=dE/hbar
        self.E_b = 0 #[eV] band gap
#        self.r = r #[m]
        self.lattice = np.array(lattice)
        self.basis = np.array(basis)
#        self.density=density #[g/cm^3]
#        self.n=self.valance*self.density*10**3/self.A #[1/m^3]
        self.n = self.valance * self.basis.shape[0] / (self.a[0]*self.a[1]*self.a[2]) * self.dz
        self.me = me #effective mass
#        self.omegaP=sp.sqrt(self.n*e**2/(eps0*m0)) #[1/s]
        self.E_p0 = hbar * sp.sqrt(self.n * e**2 / (eps0 * m0 * self.me)) #[eV]
#        self.Ep_inf = sp.sqrt(self.n*e**2/(eps0*m0))*hbar #[eV]
        self.k_F = (3 * self.n * np.pi**2)**(1/3) #[1/m]
        self.v_F = (hbar*e) *self.k_F/m0 #*e makes [m/s]        
        self.E_F = hbar**2/(2*m0) * self.k_F**2 *e #*e makes [eV]
#        self.set_Ep()
#        self.set_eps()
    
            
    def set_Ep(self, type=None, q=None):
        print('Setting E_p')
        if type is None:
            print('No type specified. Defaulting to Drude.\n',
            'Allowable values for type are:\n',
            ' -Drude: Ep(0)\n',
            ' -Lindhard: Ep(q)\n',
            ' -Lindhard2: Ep(q) with dampening')
            type = 'Drude'
            
        if type != 'Drude' and q is None:
            raise Exception('Ep can not be set when q is None, unless type=\'Drude\'.')
            
        if type=='Drude':
            """ Simple drude model at q=0"""
            self.E_p = self.E_p0
        elif type=='Lindhard':
            """Lindhard model with dampening aproaching zero"""
            q = q[:,None]
            alpha = 3/5 * self.E_F/self.Ep_0
            self.E_p = self.E_p0 + alpha* (hbar**2/m0) * q**2 *e #*e makes [eV]
        elif type=='Lindhard2':
            """Lindhard model with dampening"""
            q = q[:,None]
            alpha = 3/5 * self.E_F/self.Ep_0 * (1 - (self.E_p0/(4*self.E_F))**2)
            self.E_p = self.E_p0 + alpha* (hbar**2/m0) * q**2 *e #*e makes [eV]
            
    def set_eps(self, type=None, E=None, q=None):
        print('Setting eps')
        if E is None:
            raise Exception('eps can not be set when E is None.')
        if type is None:
            print('No type specified. Defaulting to Drude.\n',
            'Allowable values for type are:\n',
            ' -Drude: Ep(0)\n')
            type = 'Drude'
            
        if type != 'Drude' and q is None:
            raise Exception('The variable can not be set when q is None, unless type=\'Drude\'.')
        if type=='Drude':
            """ Simple drude model at q=0"""
            self.eps = 1 - self.E_p**2/(E**2-(self.E_b)**2+1j*E*self.dE)#unitless

Al=material(name="Al", a=4.046E-10, basis=[[0,0,0],[0.5, 0.5, 0],[0, 0.5, 0.5],[0.5, 0, 0.5]],
            valance=3, dE=0.66, me=1)

#da=0.1E-10#[m]
#GB=material(name="GB", a=4.046E-10+da, basis=[[0,0,0],[0.5, 0.5, 0],[0, 0.5, 0.5],[0.5, 0, 0.5]],
#            valance=3, dE=0.66, me=1)
#GB.name="GB"

def GB(da=0.1E-10,dE=0.8):
    GB=material(name="GB", a=4.046E-10+da, basis=[[0,0,0],[0.5, 0.5, 0],[0, 0.5, 0.5],[0.5, 0, 0.5]],
            valance=3, dE=0.66, me=1)
    GB.da=da
    return GB

vac=material(name="Vacuum", a=1E-9, valance=0, dE=0)

#Al2O3=material(name='Al2O3', a1=0.4046E-9, dE=10, energyB=7.6)
#
#MgO=material(name='MgO', a1=0.4046E-9, dE=6.7, energyB=7.8)
#
#Mg2Si=material(name="Mg2Si", a1=0.635E-9, valance=(2*2+4), energyB=4.6)
#
#SiO2=material(name="SiO2", a1=0.4913E-9, a2=0.5405E-9, valance=(4+2*6), energyB=8.9)
#
#Si=material(name='Si', a1=0.54307E-9, valance=2*4, dE=3.7, xtalType='Diamond')

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    E = np.arange(0.01,20,0.25)
    Al.set_Ep()
    Al.set_eps(E=E)
    vac.set_Ep()
    vac.set_eps(E=E)
    print(np.imag(vac.eps))
    
    fig, [plt1, plt2] = plt.subplots(2, 1)
    plt1.plot(E,np.imag(-1/Al.eps))
    plt2.plot(E,np.imag(-2/(Al.eps+vac.eps))+np.imag(1/Al.eps))
    plt.show()