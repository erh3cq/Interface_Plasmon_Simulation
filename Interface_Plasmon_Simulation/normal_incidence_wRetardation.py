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
m0c2=511.#[keV]
eps0=8.85E-12#[C^2/N m]

class double_differential_cross_section_normalIncidence():
    def __init__(self,microscope,spectrum,materials, t):
        #set material parameters
        self.eps = []
        for i, material in enumerate(materials):
            material.set_Ep(q=spectrum.q_perp)
            material.set_eps(E=spectrum.E)
            self.eps.append(np.nan_to_num(np.conjugate(material.eps)))
        self.t = t
        print(self.eps[0])
        self.mu2 = 1 - self.eps[0] * (microscope.v/c)**2
        
        ###Angular terms
        self.theta2 = (spectrum.q_perp/microscope.k0)**2
        self.thetaE2 = (spectrum.E/(hbar*microscope.v))**2
        self.lambda2 = self.theta2 - self.eps[0] * self.thetaE2 * (microscope.v/c)**2
        self.lambda02 = self.theta2 - self.eps[1] * self.thetaE2 * (microscope.v/c)**2
        self.phi2 = self.lambda2 + self.thetaE2
        self.phi02 = self.lambda02 + self.thetaE2
        self.phi2_01 = self.theta2 + self.thetaE2 * (1 - (self.eps[0] - self.eps[1]) * (microscope.v/c)**2)
        
        self.de = np.pi * self.t * spectrum.E / (hbar * microscope.v)
        self.Lp = np.sqrt(self.lambda02) * self.eps[0] + np.sqrt(self.lambda2) * self.eps[1] * np.tanh(np.sqrt(self.lambda2/self.thetaE2) * self.de)
        self.Lm = np.sqrt(self.lambda02) * self.eps[0] + np.sqrt(self.lambda2) * self.eps[1] / np.tanh(np.sqrt(self.lambda2/self.thetaE2) * self.de)
        
        
        self.A = 1/(np.pi * a0 * m0c2 * (microscope.v/c)**2) #[m^2/eV]
        self.B = -2 * self.theta2 * (self.eps[0] - self.eps[1])**2 / (microscope.k0 * self.phi02 * self.phi2)
        
    def bulk_mode(self):
        return self.A * np.imag(self.t * self.mu2 / (self.eps[0] * self.phi2))
    
    def surface_mode(self):
        return self.A * np.imag(self.B * self.phi2_01**4/(self.eps[0]*self.eps[1]) * 
                                                          np.sin(self.de)**2/self.Lp + np.cos(self.de)**2/self.Lm)
    def total(self):
        self.bulk_mode() + self.surface_mode()