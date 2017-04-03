# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 20:13:00 2017

@author: Eric Hoglund
"""

from __future__ import division, print_function
import numpy as np
from scipy import sqrt, imag
import scipy.integrate as integrate

import time
start_time = time.clock()

hbar=6.582E-16#[eV s]
e=1.602E-19 #[C]
hbar_J=hbar/e
c=3E8#[m/s]
m0=9.11E-31#[kg]
m0c2=511.#[keV]
eps0=8.85E-12#[C^2/N m]


class microscope():
    def __init__(self,keV=100, resolution=0.05, collection_angle=2E-3):
        self.keV=keV
        self.v=sqrt(keV*(keV+2*m0c2))/(keV+m0c2)*c #[m/s]
        self.resolution = resolution #Nion=0.05, Titan=0.05
        self.collection_angle = collection_angle #[mrad]
        self.gamma = 1/sqrt(1-self.v**2/c**2) #[au]
        self.T = 1/2*m0c2*(self.v/c)**2 #[keV]
        self.k0 = m0*self.v*self.gamma/(hbar*e) #[1/m]
    def print_parameters(self):
        print('q_beta: %.2f [1/nm]'%(self.k0*self.collection_angle/10**9))


class slab():
    def __init__(self, material, width):
        self.material = material
        self.width = width
        self.q_c = material.Ep_0 / (hbar*material.vFermi)#~1E10
        if width is not np.inf:
            self.q_min = 2*np.pi / width
        else:
            self.q_min = 0
        
        
class slab_system():
    def __init__(self, slabs, z, energy ,q_perpendicular, interface_angle=0):
        self.slabs = slabs
        self.interface_angle = interface_angle
        self.set_slab_positions()
        self.set_Ep_and_eps(energy ,q_perpendicular)
    def set_slab_positions(self):
        for index, slab in enumerate(self.slabs):
            if index==0:
                slab.position = -np.inf
            elif index==1:
                slab.position = 0
            else:
                slab.position = slabs[index-1].position + slabs[index-1].width
    def set_Ep_and_eps(self, energy ,q_perpendicular):
        for slab in self.slabs:
            slab.material.set_Ep(q=q_perpendicular)
            slab.material.set_eps(E=energy)
    def print_details(self, slab='All'):
        def _info(index):
            print('_________________')
            print('Slab ',index,':')
            print('position: %.2f [nm]   , width: %.2f [nm]'%(self.slabs[index].position/10**-9, self.slabs[index].width/10**-9))
            print('Ep_0:   %.2f [eV]   , dE:    %.2f [eV]'%(self.slabs[index].material.Ep_0, self.slabs[index].material.dE))
            print('q_min:    %.2f [1/nm] , q_c:   %.2f [1/nm]'%(self.slabs[index].q_min/10**9, self.slabs[index].q_c/10**9))    
        if slab == 'All':
            for index, slab in enumerate(self.slabs):
                _info(index)
        else:
            _info(slab)


class linescan():
    """
    energy: [eV]
    q__perpendicular: [1/m]
    xb: [m]
    """
    def __init__(self, z, energy, q_perpendicular, xb, slab_system, microscope):#introduce k_perp or theta
        self.energy = energy
        self.omega = energy / hbar
        self.z = z[:,None,None,None]
        self.xb = xb[:,None,None]
        self.slab_system = slab_system
        self.microscope = microscope
        self.q_parallel = self.omega / microscope.v
        self.q_perpendicular = np.where(q_perpendicular<=self.microscope.k0*self.microscope.collection_angle, q_perpendicular, 0)
        #self.q_perpendicular = q_perpendicular[:,None]#np.arange(k0 * theta)#Egerton pg131
        self.q = sqrt(self.q_parallel**2 + self.q_perpendicular**2)
        for index, slab in enumerate(self.slab_system):
            slab.alpha = sqrt(self.q**2 - slab.material.eps * self.omega**2/c**2)
            slab.q_perpendicular = self.q_perpendicular
            if slab.width == np.inf:
                slab.f = 1
            else:
                slab.f = np.exp(slab.alpha * slab.width)
            if slab.q_c < self.microscope.k0*self.microscope.collection_angle:
                slab.q_perpendicular=np.where(slab.q_perpendicular < slab.q_c, slab.q_perpendicular, 0)

        for index, slab in enumerate(self.slab_system):
            bulk, interface = self.slab_dP_dzdkydw(index)
            
            for key, value in enumerate(self.slice_xb(index)):
                if key != 0: bulk = np.concatenate((bulk,bulk[0][None,:]))

            if index == 0:
                self.dP_dzdkydw_bulk = bulk
                self.dP_dzdkydw_interface = interface
            else:
                self.dP_dzdkydw_bulk = np.concatenate((self.dP_dzdkydw_bulk,bulk), axis=0)
                self.dP_dzdkydw_interface =  np.concatenate((self.dP_dzdkydw_interface,interface), axis=0)
            del (bulk,interface)
            self.dP_dzdkydw = self.dP_dzdkydw_bulk + self.dP_dzdkydw_interface


    def slice_xb(self,j):
        slab_j = self.slab_system[j]
        if j==0:
            return self.xb[self.xb<0][:,None,None]
        else:
            return self.xb[np.logical_and(self.xb>=slab_j.position,self.xb<slab_j.position+slab_j.width)][:,None,None]
            
    def h(self,j,i,pm):
        slab_j = self.slab_system[j]
        slab_i = self.slab_system[i]
        return slab_j.alpha * slab_i.material.eps + pm * slab_i.alpha * slab_j.material.eps
    def ht(self,j,i,pm):
        slab_j = self.slab_system[j]
        slab_i = self.slab_system[i]
        return slab_j.alpha + pm * slab_i.alpha
    def tau(self,j):
        f = self.slab_system[j].f
        hp = self.h(j+1, j, 1)
        hm = self.h(j+1, j, -1)
        return np.array([[hp * f**2 , hm],
                         [hm * f**2 , hp]])
    def taut(self,j):
        f = self.slab_system[j].f
        htp = self.ht(j+1, j, 1)
        htm = self.ht(j+1, j, -1)
        return np.array([[htp * f**2 , htm],
                         [htm * f**2 , htp]])
    def brackets(self, n, m):
        if n==m-1:
            ar=np.array([[1,0],
                         [0,1]])
        else:
            ar=self.tau(m)
            for j in np.arange(m+1,n+1):
                ar=np.array([[self.tau(j)[0,0]*ar[0,0]+self.tau(j)[0,1]*ar[1,0] , self.tau(j)[0,0]*ar[0,1]+self.tau(j)[0,1]*ar[1,1]],
                             [self.tau(j)[1,0]*ar[0,0]+self.tau(j)[1,1]*ar[1,0] , self.tau(j)[1,0]*ar[0,1]+self.tau(j)[1,1]*ar[1,1]]])
        return ar
    def bracketst(self, n, m):
        if n==m-1:
            ar=np.array([[1,0],
                         [0,1]])
        else:
            ar=self.taut(m)
            for j in np.arange(m+1,n+1):
                ar=np.array([[self.taut(j)[0,0]*ar[0,0]+self.taut(j)[0,1]*ar[1,0] , self.taut(j)[0,0]*ar[0,1]+self.taut(j)[0,1]*ar[1,1]],
                             [self.taut(j)[1,0]*ar[0,0]+self.taut(j)[1,1]*ar[1,0] , self.taut(j)[1,0]*ar[0,1]+self.taut(j)[1,1]*ar[1,1]]])
        return ar
    def g(self, j, pm):
        slab_j = self.slab_system[j]
        xb = self.slice_xb(j)
        if j==0:
            return np.exp(pm * slab_j.alpha * xb)
        else:
            return np.exp(pm * slab_j.alpha * (xb - self.slab_system[j].position))
    def gamma(self, n, m, pm):
        return self.brackets(n,m)[0,0]*self.g(m,-1) + pm * self.brackets(n,m)[0,1]*self.g(m,1)
    def gammat(self, n, m, pm):
        return self.bracketst(n,m)[0,0]*self.g(m,-1) + pm * self.bracketst(n,m)[0,1]*self.g(m,1)
    def zeta(self, j, i, pm):
        return self.brackets(j,i)[0,0]*self.g(j+1,1) + pm * self.brackets(j,i)[1,0]*self.g(j+1,-1)
    def zetat(self, j, i, pm):
        return self.bracketst(j,i)[0,0]*self.g(j+1,1) + pm * self.bracketst(j,i)[1,0]*self.g(j+1,-1)

    def slab_chi(self,m):
        nSlabs = self.slab_system.size-2 #-1 because indexing starts at 0, and -1 because slabs is n+1
        slab_m = self.slab_system[m]
        bulkConfinment = np.pi/2 - np.arctan(slab_m.q_min/slab_m.alpha)
        
        chi_bulk = 1/slab_m.alpha * (self.microscope.v**2/c**2 - 1/slab_m.material.eps) * bulkConfinment
        chi_interface = 1/(slab_m.alpha*slab_m.material.eps * self.q**2) * (slab_m.material.eps * slab_m.q_perpendicular**2 * self.microscope.v**2/c**2 * self.zetat(m-1,0,1)*self.gammat(nSlabs,m,-1)/self.bracketst(nSlabs,0)[0,0] 
          - slab_m.alpha**2 * self.zeta(m-1,0,-1)*self.gamma(nSlabs,m,1)/self.brackets(nSlabs,0)[0,0]
          ) - 1/slab_m.alpha * (self.microscope.v**2/c**2 - 1/slab_m.material.eps)
        return chi_bulk, chi_interface
    def slab_dP_dzdkydw(self, m):
        chi_bulk, chi_interface = self.slab_chi(m)
        A = e/(4*np.pi**2*eps0*hbar**2*self.microscope.v**2)#initially e**2 but results are [1/J] so divide by e to get [1/eV]
        
        dP_dzdkydw_bulk = A*imag(chi_bulk[None,:,:])
        dP_dzdkydw_interface = A*imag(chi_interface)
        return dP_dzdkydw_bulk, dP_dzdkydw_interface
    def slab_dP_dzdw(self):
        dP_dzdw_bulk = integrate.simps(self.dP_dzdkydw_bulk, self.q_perpendicular[None,:], axis=-2)
        dP_dzdw_interface = integrate.simps(self.dP_dzdkydw_interface, self.q_perpendicular[None,:], axis=-2)
        return dP_dzdw_bulk+dP_dzdw_interface, dP_dzdw_bulk, dP_dzdw_interface
    def slab_dP_dkydw(self):
        dP_dkydw_bulk = self.dP_dzdkydw_bulk*self.z.max()#slab_m.dP_dzdkydw_bulk, self.z, axis=-4)
        dP_dkydw_interface = self.dP_dzdkydw_interface*self.z.max()#integrate.simps(slab_m.dP_dzdkydw_interface, self.z, axis=-4)
        return dP_dkydw_bulk+dP_dkydw_interface, dP_dkydw_bulk, dP_dkydw_interface


microscope = microscope()
microscope.print_parameters()

thickness = 20E-10
z = np.linspace(0, thickness, num=50)

minE, maxE=5, 25
energy = np.arange(minE, maxE, microscope.resolution)
q_perpendicular = np.linspace(0,microscope.k0*microscope.collection_angle,200)[:,None]#np.arange(0,1E10,1E5)


from materials import Al, GB, vac
slabs = np.array([slab(Al, np.inf) , slab(GB(),0.5*4/sqrt(3) * 10**-10) , slab(Al, np.inf)])
#slabs = np.array([slab(vac, np.inf) , slab(Al,np.inf)])
#slabs = np.array([slab(Al,np.inf)])
slabs = slab_system(slabs, z, energy ,q_perpendicular)
slabs.print_details(slab='All')

inSlab = 0
scanLen = np.arange(-5,5,0.1)*10**-9#np.array([-20,-0.5,-0.09,0.09,0.5,100]) * 10**-10



lossfunction = linescan(z, energy, q_perpendicular, scanLen, slabs.slabs, microscope)

print("Calculation time: %s min ---" % ((time.clock() - start_time)/60))

import matplotlib.pyplot as plt

bounds = [minE,maxE,q_perpendicular.min() *10**-9,q_perpendicular.max() *10**-9]
total, bulk, interface = lossfunction.slab_dP_dkydw()

fig_dP_dkydw = plt.figure('dP2_dqy dE)')
plt.subplot(2,2,1)
plt.imshow(total[0] /10**-9, origin='lower', aspect='auto', extent=bounds)
plt.colorbar().set_label(label=r'$\frac{dP^2}{dq_y\ dE}\ \left[ \frac{nm}{eV} \right]$', fontsize=16)
plt.title('Total')
plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')
plt.ylabel(r'$q_\perp [nm^-]$', fontsize=12, fontweight='bold')

plt.subplot(2,2,2)
plt.imshow(bulk[0] /10**-9, origin='lower', aspect='auto', extent=bounds)
plt.colorbar().set_label(label=r'$\frac{dP^2}{dq_y\ dE}\ \left[ \frac{nm}{eV} \right]$', fontsize=16)
plt.title('Bulk')
plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')
plt.ylabel(r'$q_\perp [nm^-]$', fontsize=12, fontweight='bold')

plt.subplot(2,2,3)
plt.imshow(interface[0] /10**-9, origin='lower', aspect='auto', extent=bounds)
plt.colorbar().set_label(label=r'$\frac{dP^2}{dq_y\ dE}\ \left[ \frac{nm}{eV} \right]$', fontsize=16)
plt.title('Interface')
plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')
plt.ylabel(r'$q_\perp [nm^-]$', fontsize=12, fontweight='bold')

plt.get_current_fig_manager().window.showMaximized()



bounds = [minE,maxE,q_perpendicular.min() *10**-9,q_perpendicular.max() *10**-9]

for pos in np.array([0,round((scanLen.size-1)/2),scanLen.size-1]):
    plt.figure('dP3_dz dqy dE at %d '%(pos))
    totalplot = plt.subplot(2,2,1)
    plt.imshow(lossfunction.dP_dzdkydw[pos], origin='lower', aspect='auto', extent=bounds)
    plt.colorbar().set_label(label=r'$\frac{dP^3}{dz\ dq_y\ dE}\ [eV^-]$', fontsize=16)
    plt.title('Total')
    plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')
    plt.ylabel(r'$q_\perp [nm^-]$', fontsize=12, fontweight='bold')
    
    bulklplot = plt.subplot(2,2,2)
    plt.imshow(lossfunction.dP_dzdkydw_bulk[pos], origin='lower', aspect='auto', extent=bounds)
    plt.colorbar().set_label(label=r'$\frac{dP^3}{dz\ dq_y\ dE}\ [eV^-]$', fontsize=16)
    plt.title('Bulk')
    plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')
    plt.ylabel(r'$q_\perp [nm^-]$', fontsize=12, fontweight='bold')
    
    interfaceplot = plt.subplot(2,2,3)
    plt.imshow(lossfunction.dP_dzdkydw_interface[pos], origin='lower', aspect='auto', extent=bounds)
    plt.colorbar().set_label(label=r'$\frac{dP^3}{dz\ dq_y\ dE}\ [eV^-]$', fontsize=16)
    plt.title('Interface')
    plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')
    plt.ylabel(r'$q_\perp [nm^-]$', fontsize=12, fontweight='bold')
    
    plt.get_current_fig_manager().window.showMaximized()


bounds = [minE,maxE,scanLen.min()/10**-9,scanLen.max()/10**-9]
total, bulk, interface = lossfunction.slab_dP_dzdw()

fig_dP_dzdw=plt.figure('dP2_dz dE',)
totalplot = plt.subplot(2,2,1)
plt.imshow(total *10**-9, origin='lower', aspect='auto' ,extent=bounds)
plt.colorbar().set_label(label=r'$\frac{dP^2}{dz\ dE}\ \left[ \frac{1}{eV\ nm}\right] $', fontsize=16)
plt.title('Total')
plt.ylabel(r'$x [nm^-]$', fontsize=12, fontweight='bold')
plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')

bulklplot = plt.subplot(2,2,2)
plt.imshow(bulk *10**-9, origin='lower', aspect='auto'
                     ,extent=bounds)
plt.colorbar().set_label(label=r'$\frac{dP^2}{dz\ dE}\ \left[ \frac{1}{eV\ nm}\right] $', fontsize=16)
plt.title('Bulk')
plt.ylabel(r'$x [nm^-]$', fontsize=12, fontweight='bold')
plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')

interfaceplot = plt.subplot(2,2,3)
plt.imshow(interface *10**-9, origin='lower', aspect='auto'
                     ,extent=bounds)
plt.colorbar().set_label(label=r'$\frac{dP^2}{dz\ dE}\ \left[ \frac{1}{eV\ nm}\right] $', fontsize=16)
plt.title('Interface')
plt.ylabel(r'$x [nm^-]$', fontsize=12, fontweight='bold')
plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')

plt.get_current_fig_manager().window.showMaximized()




plt.show()