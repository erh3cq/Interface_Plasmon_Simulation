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
eps0=8.85E-12#[C^2/N m^2]
a0 = 4*np.pi*eps0 * hbar**2 / m0 #4*pi*eps0 *hbar^2/(m0*e^2) where hbar is in Js


class material():
    #TODO: add effective mass
    def __init__(self, name='General Material', a=1E-9, dz=1,
                 lattice = [[1, 0, 0],[0, 1, 0],[0, 0, 1]],
                 basis = [[0,0,0]],
                 valance=3, dE=10, me=1, eps_inf=1):
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
        self.na = self.basis.shape[0] / (self.a[0]*self.a[1]*self.a[2]) * self.dz
        self.n = self.valance * self.na
        self.me = me #effective mass
        self.eps_inf = eps_inf
#        self.omegaP=sp.sqrt(self.n*e**2/(eps0*m0)) #[1/s]
        self.E_p0 = hbar * sp.sqrt(self.n * e**2 / (eps0 * m0 * self.me)) #[eV]
#        self.Ep_inf = sp.sqrt(self.n*e**2/(eps0*m0))*hbar #[eV]
        self.k_F = (3 * self.n * np.pi**2)**(1/3) #[1/m]
        self.v_F = (hbar*e) *self.k_F/(m0*me) #*e makes [m/s]        
        self.E_F = hbar**2/(2*m0*me) * self.k_F**2 *e #*e makes [eV]
        self.q_c = (m0*self.me)/(hbar* hbar*e) * self.E_p0/self.k_F
    
            
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
            self.E_p = self.E_p0 + alpha* (hbar**2/(m0*self.me)) * q**2 *e #*e makes [eV]
        elif type=='Lindhard2':
            """Lindhard model with dampening"""
            q = q[:,None]
            alpha = 3/5 * self.E_F/self.Ep_0 * (1 - (self.E_p0/(4*self.E_F))**2)
            self.E_p = self.E_p0 + alpha* (hbar**2/(m0*self.me)) * q**2 *e #*e makes [eV]
            
    def set_eps(self, type=None, E=None, q=None, multiple_eps=False):
        print('Setting eps')
        if E is None:
            raise Exception('eps can not be set when E is None.')
        elif q is None and type=='Lindhard':
            raise Exception('eps_L can not be set when q is None.')
        if type is None:
            print('No type specified. Defaulting to Drude.\n',
            'Allowable values for type are:\n',
            ' -Drude: Ep(0)\n')
            type = 'Drude'
            
        if type != 'Drude' and q is None:
            raise Exception('The variable can not be set when q is None, unless type=\'Drude\'.')
        if type=='Drude':
            """ Simple drude model at q=0"""
            self.eps = self.eps_inf + self.E_p**2/(self.E_b**2 - E**2 - 1j*E*self.dE)#unitless
            if multiple_eps is True:
                self.eps_D = self.eps
        
        if type=='Lindhard':
            """ Lindhard model from Abajo (2010)"""
            def Rxy(val1,val2):
#                print('test1',1/(2*val1))
#                print('test2',(1-((val1**2+val2)/(2*val2))**2))
#                print('test3',np.log((val1**2+2*val1+val2)/(val1**2-2*val1+val2)))
                return np.nan_to_num(1/(2*val1)) * (1-((val1**2+val2)/(2*val1))**2) * np.log((val1**2+2*val1+val2)/(val1**2-2*val1+val2))
            E = (E + 1j*self.dE)
            self.eps = self.eps_inf + (m0*self.me)/hbar**2 * 1/(2*eps0*np.pi**2) * self.k_F/q**2 *\
            (1 + Rxy(q/self.k_F, E/self.E_F) + Rxy(q/self.k_F, -1*E/self.E_F))#unitless
            if multiple_eps is True:
                self.eps_L = self.eps
        
        if type=='Mermin':
            """ Mermin model"""
            Ei = (E + 1j*self.dE)
            def Rxy(val1,val2):
#                print('test1',1/(2*val1))
#                print('test2',(1-((val1**2+val2)/(2*val2))**2))
#                print('test3',np.log((val1**2+2*val1+val2)/(val1**2-2*val1+val2)))
                return np.nan_to_num(1/(2*val1)) * (1-((val1**2+val2)/(2*val1))**2) * np.log((val1**2+2*val1+val2)/(val1**2-2*val1+val2))
            def Chi_L():
                return (m0*self.me)/hbar**2 * 1/(2*eps0*np.pi**2) * self.k_F/q**2 *\
                (1 + Rxy(q/self.k_F, Ei/self.E_F) + Rxy(q/self.k_F, -1*Ei/self.E_F))#unitless
            def Chi_L0():
                x = q/self.k_F
                return (m0*self.me)/hbar**2 * 1/(2*eps0*np.pi**2) * self.k_F/q**2 *\
                (1 + (4-x**2)/(4*x) * np.log(np.abs((x+2)/(x-2))))#unitless
            self.eps = self.eps_inf + Ei*Chi_L() / (E + 1j*self.dE * Chi_L()/Chi_L0())            
            if multiple_eps is True:
                self.eps_M = self.eps

Al=material(name="Al", a=4.046E-10, basis=[[0,0,0],[0.5, 0.5, 0],[0, 0.5, 0.5],[0.5, 0, 0.5]],
            valance=3, dE=0.66, me=1.025, eps_inf=1.037)

#da=0.1E-10#[m]
#GB=material(name="GB", a=4.046E-10+da, basis=[[0,0,0],[0.5, 0.5, 0],[0, 0.5, 0.5],[0.5, 0, 0.5]],
#            valance=3, dE=0.66, me=1)
#GB.name="GB"

def GB(da=0.1E-10,dE=0.8):
    GB=material(name="GB", a=4.046E-10+da, basis=[[0,0,0],[0.5, 0.5, 0],[0, 0.5, 0.5],[0.5, 0, 0.5]],
            valance=3, dE=dE, me=1.025, eps_inf=1.037)
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
    np.seterr(divide='ignore')
    
    q = np.linspace(0,1.5,200,endpoint=True)*1E10
    E = np.arange(13,17,0.03)

    q = q[:,np.newaxis]
    E = E[np.newaxis,:]

    materials = [Al, GB()]
    
    for material in materials:
        print(material.name)
        print('n: ',material.n,'[e-/m^3]')
        print('k_F: ',material.k_F)
        print('E_F: ',material.E_F)
        material.set_Ep(type='Drude')
        print('Ep_0: ',material.E_p0)
        print('q_c:',material.q_c)
        print(material.E_p0 + 3/5 * material.E_F/material.E_p0 * (hbar**2/(m0*material.me)) * (material.q_c)**2 *e)
        print(3/5 * material.E_F/material.E_p0)
        #material.set_eps(E=E, q=q, type='Lindhard', multiple_eps=True)
        #print('eps:')
        #print(material.eps_L)
        print('')
#        print('test',np.nan_to_num(np.inf)*np.array([0,1,2]))
        
    #fig, [plt1, plt2] = plt.subplots(2, 1)
    #plt1.plot(E,np.imag(-1/Al.eps))
    #plt2.plot(E,np.imag(-2/(Al.eps+vac.eps))+np.imag(1/Al.eps))
    #plt.show()