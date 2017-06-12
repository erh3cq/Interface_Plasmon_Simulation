# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 19:06:40 2017

@author: Eric Hoglund


Egerton pg 162

Generate an Energy loss spectra with
-Volume Plasmon (Cˇerenkov-enhanced)
-Surface Plasmon
-Guided-light modes

"""

from __future__ import division, print_function
import numpy as np

hbar=6.582E-16#[eV s]
e=1.602E-19 #[C]
a0 = 5.292*10**-11 #[m]
c=3E8#[m/s]
m0=9.11E-31#[kg]
m0c2=511E3#[eV]
eps0=8.85E-12#[C^2/N m]

class double_differential_cross_section_normalIncidence():
    def __init__(self, microscope, spectrum, materials, t):
        self.beta2 = (microscope.v/c)**2
        self.q_parallel = spectrum.E/(hbar*microscope.v)
        
        ###set material parameters
        self.eps = []
        for i, material in enumerate(materials):
            material.set_Ep(q=np.sqrt(spectrum.q_perpendicular**2 + self.q_parallel**2))
            material.set_eps(E=spectrum.E)
            self.eps.append(np.conjugate(material.eps))
            print(material.name)
            print('Ep_0',material.Ep_0,'[eV]')
            print('E_fermi',material.E_fermi,'[eV]')
            print()
        self.t = t
        
        self.mu2 = 1 - self.eps[1] * self.beta2
        
        ###Angular terms
        self.theta2 = (spectrum.q_perpendicular/microscope.k0)**2
        self.thetaE2 = (spectrum.E/(2*microscope.gamma*microscope.T))**2 #(spectrum.E/(hbar*microscope.v*microscope.k0))**2

        self.lambda2 = self.theta2 - self.eps[1] * self.thetaE2 * self.beta2
        self.lambda02 = self.theta2 - self.eps[0] * self.thetaE2 * self.beta2
        self.phi2 = self.lambda2 + self.thetaE2
        self.phi02 = self.lambda02 + self.thetaE2
        self.phi2_01 = self.theta2 + self.thetaE2 * (1 - (self.eps[1] + self.eps[0]) * self.beta2)
        
        self.de = self.t * spectrum.E / (2*hbar * microscope.v)
        self.tanh = np.tanh(np.sqrt(self.lambda2/self.thetaE2) * self.de)
        self.Lp = np.sqrt(self.lambda02) * self.eps[1] + np.sqrt(self.lambda2) * self.eps[0] * self.tanh #symetric mode     
        self.Lm = np.sqrt(self.lambda02) * self.eps[1] + np.sqrt(self.lambda2) * self.eps[0] / self.tanh #anti-symetric mode
        
        
        self.A =e/(np.pi**2 * a0 * m0 * microscope.v**2)  #e**3 / (4*np.pi**3 * hbar**2 * eps0 * microscope.v**2) #[m^2/eV]
        self.B = -2 * self.theta2 * (self.eps[1] - self.eps[0])**2 / (microscope.k0 * (self.phi02 * self.phi2)**2)
        
    def bulk_mode(self):
        return self.A * np.imag(self.t * self.mu2 / (self.eps[1] * self.phi2))
    
    def surface_mode(self):
        return self.A * np.imag(self.B * self.phi2_01**2 / (self.eps[1]*self.eps[0])
                               * (np.sin(self.de)**2/self.Lp + np.cos(self.de)**2/self.Lm))
    def guidedLight1(self):
        return self.A * np.imag(self.B * self.beta2 * np.sqrt(self.lambda02*self.thetaE2) * self.phi2_01 / self.eps[0] *
                                (1/self.Lp - 1/self.Lm) * np.sin(2 * self.de))
        
    def guidedLight2(self):
        return self.A * np.imag(self.B * -1* self.beta2**2 * np.sqrt(self.lambda02*self.lambda2) * self.thetaE2 *
                                (np.cos(self.de)**2 * self.tanh/self.Lp + np.sin(self.de)**2/(self.Lm * self.tanh)))
        
    def total(self):
        return self.bulk_mode() + self.surface_mode() + self.guidedLight1() + self.guidedLight2()