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
        self.q_c = material.Ep_inf / (hbar*material.vFermi)#~1E10
        if width is not np.inf:
            self.q_min = 2*np.pi / width
        else:
            self.q_min = 0
        
        
class slab_system():
    def __init__(self, slabs, z, energy ,q_perpendicular, interface_angle=2E-3):
        self.slabs = slabs
        self.interface_angle = interface_angle
        self.set_slab_positions()
        self.set_Ep_and_eps(energy ,q_perpendicular)
    def set_slab_positions(self):
        for index, slab in enumerate(self.slabs):
            if index==0:
                slab.position = -np.inf - z * self.interface_angle
            elif index==1:
                slab.position = 0 - z * self.interface_angle
            else:
                slab.position = slabs[index-1].position + slabs[index-1].width# - z * self.interface_angle
    def set_Ep_and_eps(self, energy ,q_perpendicular):
        for slab in self.slabs:
            slab.material.set_Ep(q=q_perpendicular)
            slab.material.set_eps(E=energy)
    def print_details(self, slab='All'):
        def _info(index):
            print('_________________')
            print('Slab ',index,':')
            print('position: %.2f [nm]   , width: %.2f [nm]'%(self.slabs[index].position[0]/10**-9, self.slabs[index].width/10**-9))
            print('Ep_inf:   %.2f [eV]   , dE:    %.2f [eV]'%(self.slabs[index].material.Ep_inf, self.slabs[index].material.dE))
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
        self.z = z
        self.xb = xb[:, None, None]
        self.slab_system = slab_system
        self.microscope = microscope
        self.q_parallel = self.omega / microscope.v
        self.q_perpendicular = np.where(q_perpendicular<=self.microscope.k0*self.microscope.collection_angle, q_perpendicular, 0)[:,None]
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
        
#        self.dP_dzdkydw = np.empty((self.z.size, self.xb.size, self.q_perpendicular.size, self.energy.size))
#        self.dP_dzdkydw_bulk = np.empty_like(self.dP_dzdkydw)
#        self.dP_dzdkydw_interface = np.empty_like(self.dP_dzdkydw)
        for z_index, z_value in enumerate(self.z):
            self.z_pos = z_index
            print('z=',z_index)
            for index, slab in enumerate(self.slab_system):
                print('slab:',index)
#                self.set_slab_chi(index)
                bulk, interface = self.slab_dP_dzdkydw(index)
                
                for key, value in enumerate(self.slice_xb(index, z_pos=0)):
                    if key != 0: bulk = np.concatenate((bulk,bulk[0][None,:]))
                    
#                dP_dzdkydw_bulk_temp = np.empty_like(slab.dP_dzdkydw_interface)
#                for key, value in enumerate(self.slice_xb(index, z_pos=0)):
#                    dP_dzdkydw_bulk_temp[key, :,:] =  slab.dP_dzdkydw_bulk
    
                if index == 0:
#                    dP_dzdkydw_bulk_temp2 = dP_dzdkydw_bulk_temp[None, :]
                    dP_dzdkydw_bulk_temp = bulk[None, :]
                    print('b shape:',dP_dzdkydw_bulk_temp.shape)
#                    self.dP_dzdkydw_bulk = bulk[None, :]
#                    dP_dzdkydw_interface_temp = slab.dP_dzdkydw_interface[None, :]
                    dP_dzdkydw_interface_temp = interface[None, :]
                    print('i shape:',dP_dzdkydw_interface_temp.shape)
#                    self.dP_dzdkydw_interface = interface[None, :]
                else:
#                    dP_dzdkydw_bulk_temp2 = np.concatenate((dP_dzdkydw_bulk_temp2,dP_dzdkydw_bulk_temp[None, :]), axis=1)
                    dP_dzdkydw_bulk_temp = np.concatenate((dP_dzdkydw_bulk_temp,bulk[None, :]), axis=1)
                    print('b shape:',dP_dzdkydw_bulk_temp.shape)
#                    self.dP_dzdkydw_bulk = np.concatenate((self.dP_dzdkydw_bulk,bulk[None, :]), axis=1)
#                    dP_dzdkydw_interface_temp =  np.concatenate((dP_dzdkydw_interface_temp,slab.dP_dzdkydw_interface[None, :]), axis=1)
                    dP_dzdkydw_interface_temp =  np.concatenate((dP_dzdkydw_interface_temp,interface[None, :]), axis=1)
                    print('i shape:',dP_dzdkydw_interface_temp.shape)
#                    self.dP_dzdkydw_interface =  np.concatenate((self.dP_dzdkydw_interface,interface[None, :]), axis=1)
            if z_index == 0:
                self.dP_dzdkydw_bulk = dP_dzdkydw_bulk_temp
                print('b',self.dP_dzdkydw_bulk.shape)
                self.dP_dzdkydw_interface = dP_dzdkydw_interface_temp
                print('i',self.dP_dzdkydw_interface.shape)
            else:
                self.dP_dzdkydw_bulk = np.concatenate((self.dP_dzdkydw_bulk,dP_dzdkydw_bulk_temp), axis=0)
                print('b',self.dP_dzdkydw_bulk.shape)
                self.dP_dzdkydw_interface = np.concatenate((self.dP_dzdkydw_interface,dP_dzdkydw_interface_temp), axis=0)
                print('i',self.dP_dzdkydw_interface.shape)
            del dP_dzdkydw_bulk_temp
#            del dP_dzdkydw_bulk_temp2
            del dP_dzdkydw_interface_temp
            self.dP_dzdkydw = self.dP_dzdkydw_bulk + self.dP_dzdkydw_interface
            
                
            
    def slice_xb(self, j, z_pos=0):
        slab_j = self.slab_system[j]
#        print('xb shape ',self.xb.shape)
#        print('pos size ',slab_j.position.shape)
        if j==0:
#            print('cond ',np.where(self.xb<0,self.xb,np.nan).shape)
            #print('end shape ',self.xb[self.xb<0][:,None,None,None].shape)
            return self.xb[self.xb<0][:,None,None]
        else:
            #print('cond ',np.logical_and(self.xb>=slab_j.position,self.xb<slab_j.position+slab_j.width).shape)
            #print('end shape ',self.xb[np.logical_and(self.xb>=slab_j.position,self.xb<slab_j.position+slab_j.width)][:,None,None,None].shape)
            return self.xb[np.logical_and(self.xb>slab_j.position[z_pos],self.xb<=slab_j.position[z_pos]+slab_j.width)][:,None,None]
            
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
        xb = self.slice_xb(j, z_pos=0)
        if j==0:
            return np.exp(pm * slab_j.alpha * xb)
        else:
#            print('xb',(xb).shape)
#            print('pos',(xb - self.slab_system[j].position).shape)
#            print('-',(self.slab_system[j].position).shape)
            return np.exp(pm * slab_j.alpha * (xb - self.slab_system[j].position[0]))
    def gamma(self, n, m, pm):
        return self.brackets(n,m)[0,0]*self.g(m,-1) + pm * self.brackets(n,m)[0,1]*self.g(m,1)
    def gammat(self, n, m, pm):
        return self.bracketst(n,m)[0,0]*self.g(m,-1) + pm * self.bracketst(n,m)[0,1]*self.g(m,1)
    def zeta(self, j, i, pm):
        return self.brackets(j,i)[0,0]*self.g(j+1,1) + pm * self.brackets(j,i)[1,0]*self.g(j+1,-1)
    def zetat(self, j, i, pm):
        return self.bracketst(j,i)[0,0]*self.g(j+1,1) + pm * self.bracketst(j,i)[1,0]*self.g(j+1,-1)
    
    def slab_chi(self, m):
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
        
    def set_slab_dP_dzdw(self, m):
        slab_m = self.slab_system[m]
        slab_m.dP_dzdw_bulk = integrate.simps(slab_m.dP_dzdkydw_bulk, slab_m.q_perpendicular, axis=-2)
        slab_m.dP_dzdw_interface = integrate.simps(slab_m.dP_dzdkydw_interface, slab_m.q_perpendicular[None,None,:], axis=-2)
        slab_m.dP_dzdw = slab_m.dP_dzdw_bulk + slab_m.dP_dzdw_interface
    def set_slab_dP_dkydw(self, m):
        slab_m = self.slab_system[m]
        slab_m.dP_dkydw_bulk = slab_m.dP_dzdkydw_bulk*self.z.max()#slab_m.dP_dzdkydw_bulk, self.z, axis=-4)
        slab_m.dP_dkydw_interface = slab_m.dP_dzdkydw_interface*self.z.max()#integrate.simps(slab_m.dP_dzdkydw_interface, self.z, axis=-4)
        slab_m.dP_dkydw = slab_m.dP_dkydw_bulk + slab_m.dP_dkydw_interface
        
    def slab_dP_dzdw(self):
        print(self.dP_dzdkydw_bulk.shape)
        dP_dzdw_bulk = integrate.simps(self.dP_dzdkydw_bulk, self.q_perpendicular[None,None,:], axis=-2)
        dP_dzdw_interface = integrate.simps(self.dP_dzdkydw_interface, self.q_perpendicular[None,None,:], axis=-2)
        dP_dzdw = dP_dzdw_bulk + dP_dzdw_interface
        return dP_dzdw, dP_dzdw_bulk, dP_dzdw_interface
    def slab_dP_dkydw(self):
        print(self.dP_dzdkydw_bulk.shape, self.z.shape)
        dP_dkydw_bulk = integrate.simps(self.dP_dzdkydw_bulk, self.z, axis=0)#slab_m.dP_dzdkydw_bulk, self.z, axis=-4)
        dP_dkydw_interface = integrate.simps(self.dP_dzdkydw_interface, self.z, axis=0) #integrate.simps(slab_m.dP_dzdkydw_interface, self.z, axis=-4)
        dP_dkydw = dP_dkydw_bulk + dP_dkydw_interface
        return dP_dkydw, dP_dkydw_bulk, dP_dkydw_interface


microscope = microscope()
microscope.print_parameters()

thickness = 20E-10
z = np.linspace(0, thickness, num=6)[:, None, None, None]

minE, maxE=5, 25
energy = np.arange(minE, maxE, microscope.resolution)
q_perpendicular = np.linspace(0,microscope.k0*microscope.collection_angle,200)#np.arange(0,1E10,1E5)


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

for pos in np.array([0,round((scanLen.size-1)/2),scanLen.size-1]):
    plt.figure('dP^3/(dz dq_y dE) at %d '%(pos))
    totalplot = plt.subplot(2,2,1)
    plt.imshow(lossfunction.dP_dzdkydw[0, pos], origin='lower', aspect='auto', extent=bounds)
    plt.colorbar(label=r'$\frac{dP^3}{dz dq_y dE} [eV^-]$')
    plt.title('Total')
    plt.xlabel(r'$E [eV]$')
    plt.ylabel(r'$q_\perp [nm^-]$')
    
    bulklplot = plt.subplot(2,2,2)
    plt.imshow(lossfunction.dP_dzdkydw_bulk[0, pos], origin='lower', aspect='auto', extent=bounds)
    plt.colorbar(label=r'$\frac{dP^3}{dz dq_y dE} [eV^-]$')
    plt.title('Bulk')
    plt.xlabel(r'$E [eV]$')
    plt.ylabel(r'$q_\perp [nm^-]$')
    
    interfaceplot = plt.subplot(2,2,3)
    plt.imshow(lossfunction.dP_dzdkydw_interface[0, pos], origin='lower', aspect='auto', extent=bounds)
    plt.colorbar(label=r'$\frac{dP^3}{dz dq_y dE} [eV^-]$')
    plt.title('Interface')
    plt.xlabel(r'$E [eV]$')
    plt.ylabel(r'$q_\perp [nm^-]$')


    
total, bulk, interface = lossfunction.slab_dP_dkydw()
bounds = [minE,maxE,q_perpendicular.min() *10**-9,q_perpendicular.max() *10**-9]
fig_dP_dkydw = plt.figure('dP^3/(dq_y dE) at %d '%(pos),)
plt.subplot(2,2,1)
plt.imshow(total[0] /10**-9, origin='lower', aspect='auto', extent=bounds)
plt.colorbar().set_label(label=r'$\frac{dP^2}{dq_y\ dE}\ \left[ \frac{nm}{eV} \right]$', fontsize=16)
plt.title('Total')
plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')
plt.ylabel(r'$q_\perp [nm^-]$', fontsize=12, fontweight='bold')

plt.subplot(2,2,2)
plt.imshow(bulk[0] /10**-9, origin='lower', aspect='auto', extent=bounds)
plt.colorbar().set_label(label=r'$\frac{dP^2}{dq_y\ dE}\ \left[ \frac{nm}{eV} \right]$', fontsize=16)
plt.title('Total')
plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')
plt.ylabel(r'$q_\perp [nm^-]$', fontsize=12, fontweight='bold')

plt.subplot(2,2,3)
plt.imshow(interface[0] /10**-9, origin='lower', aspect='auto', extent=bounds)
plt.colorbar().set_label(label=r'$\frac{dP^2}{dq_y\ dE}\ \left[ \frac{nm}{eV} \right]$', fontsize=16)
plt.title('Total')
plt.xlabel(r'$E [eV]$', fontsize=12, fontweight='bold')
plt.ylabel(r'$q_\perp [nm^-]$', fontsize=12, fontweight='bold')


total, bulk, interface = lossfunction.slab_dP_dzdw()
bounds = [minE,maxE,scanLen.min()/10**-9,scanLen.max()/10**-9]    
fig_dP_dzdw=plt.figure('dP^3/(dz dE) at z=0')
totalplot = plt.subplot(2,2,1)
plt.imshow(total[0] *10**-9, origin='lower', aspect='auto' ,extent=bounds)
plt.colorbar(label=r'$\frac{dP^3}{dz dE} \left[ \frac{1}{eV nm}\right] $')
plt.title('Total')
plt.ylabel(r'$x [nm^-]$')
plt.xlabel(r'$E [eV]$')

bulklplot = plt.subplot(2,2,2)
plt.imshow(bulk[0] *10**-9, origin='lower', aspect='auto'
                     ,extent=bounds)
plt.colorbar(label=r'$\frac{dP^3}{dz dE} \left[ \frac{1}{eV nm}\right] $')
plt.title('Bulk')
plt.ylabel(r'$x [nm^-]$')
plt.xlabel(r'$E [eV]$')

interfaceplot = plt.subplot(2,2,3)
plt.imshow(interface[0] *10**-9, origin='lower', aspect='auto'
                     ,extent=bounds)
plt.colorbar(label=r'$\frac{dP^3}{dz dE} \left[ \frac{1}{eV nm}\right] $')
plt.title('Interface')
plt.ylabel(r'$x [nm^-]$')
plt.xlabel(r'$E [eV]$')
plt.show()