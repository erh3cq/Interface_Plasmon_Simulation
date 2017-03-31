# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 17:11:43 2016

@author: Eric Hoglund
"""

from __future__ import division
import numpy as np
import scipy as sp

hbar=6.582E-16#[eV s]   1.055E-34[J s]
u=1.660539040E-27
e=1.6E-19 #[C]
m0=9.11E-31#[kg]
m0c2=511E3#[eV]
eps0=8.85E-12#[C^2/N m]

def xtal(xtalType):
    if xtalType=='FCC':
        lattice = np.array([[0,1/2,1/2],[1/2,1/2,0],[1/2,0,1/2]])
        basis = np.array([[0,0,0]])
#        basis=np.array([[0,0,0],[0,1/2,1/2],[1/2,1/2,0],[1/2,0,1/2]])
#        lattice=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
    elif xtalType=='Diamond':
        basis=np.array([[0,0,0],[1/4,1/4,1/4]])
        lattice=np.array([[0,0,0],[0,1/2,1/2],[1/2,1/2,0],[1/2,0,1/2]])
    return lattice, basis;


class material():
    def __init__(self,name='General Material', a1=1E-9, a2=None, a3=None, dz=1, valance=3, dE=10, energyB=0, r=0.118E-9, xtalType='FCC', density=2.702, mass=26.9815385):
        self.name=name
        self.a1=a1#[m]
        if a2 is None: self.a2=a1
        else: self.a2=a2
        if a3 is None: self.a3=a1
        else: self.a3=a3
        self.a=np.array([self.a1,self.a2,self.a3])
        self.dz=dz #density adjustment factor (ex. on average 1 of 4 atoms in FCC is missing then dz=0.75)
        self.valance=valance #valance per basis
        self.dE=dE
        self.damp=dE/hbar
        self.energyB=energyB #[eV] band gap
        self.r=r #[m]
        self.lattice, self.basis=xtal(xtalType)
        self.density=density #[g/cm^3]
        self.A=mass*u#[kg]
        #self.n=self.valance*self.density*10**3/self.A #[1/m^3]
        self.na = (np.size(self.lattice,axis=1)+1)/(self.a1*self.a2*self.a3) #atoms/[m]^3
        self.n=self.valance*self.na*self.dz #e-/[m]^3
        self.omegaP=sp.sqrt(self.n*e**2/(eps0*m0)) #[1/s]
        self.energyP=sp.sqrt(self.n*e**2/(eps0*m0))*hbar #[eV]
        self.Ep_0 = sp.sqrt(self.n*e**2/(eps0*m0))*hbar #[eV]
        self.kFermi=(3*self.n*np.pi**2)**(1/3) #[1/m]
        self.E_fermi=hbar**2/(2*m0) * self.kFermi**2 *e #*e makes [eV]
        self.vFermi=hbar*self.kFermi/m0 *e #*e makes [m/s]
        self.set_Ep()
        self.set_eps()
    def set_Ep(self, q=None):
        if q is None:
            self.Ep = self.Ep_0
        else:
            q = q#[:,None]
            alpha = np.nan_to_num(3* self.E_fermi / (5*self.Ep_0) * (1-(self.Ep_0/(4*self.E_fermi))**2))
            #alpha = np.nan_to_num(3* self.E_fermi / (5*self.Ep_0))
            self.Ep = self.Ep_0 + alpha* ((hbar)**2*e/m0) * q**2 #*e makes [eV]
            print('A',(hbar)**2*e/m0)
            print('alpha',alpha)
    def set_eps(self, E=None):
        if E is not None:
            self.eps = 1 - self.Ep**2/(E**2-(self.energyB)**2+1j*E*self.damp*hbar)#unitless

Al=material(name="Al", a1=0.4046E-9, valance=3, dE=0.66, r=0.118E-9, density=2.702, mass=26.9815385)

#da=0.1E-10#[m]
#GB=material(Al, a1=Al.a1+da, a2=Al.a1, a3=Al.a3, dE=5)
#GB.name="GB"

def GB(da=0.1E-10,dE=0.8):
    GB=material(name="GB", a1=Al.a1+da, a2=Al.a1, a3=Al.a3, dE=dE)
    GB.da=da
    return GB

vac=material(name="Vacuum", a1=1E-9, valance=0, dE=0)

Al2O3=material(name='Al2O3', a1=0.4046E-9, dE=22, energyB=7.6)
Al2O3.Ep_0 = 24

MgO=material(name='MgO', a1=0.4046E-9, dE=6.7, energyB=7.8)

Mg2Si=material(name="Mg2Si", a1=0.635E-9, valance=(2*2+4), energyB=4.6)

SiO2=material(name="SiO2", a1=0.4913E-9, a2=0.5405E-9, valance=(4+2*6), energyB=8.9)

Si=material(name='Si', a1=0.54307E-9, valance=2*4, dE=3.7, xtalType='Diamond')