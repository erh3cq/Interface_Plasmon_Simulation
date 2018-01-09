# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 18:36:58 2017

@author: erhog
"""

from __future__ import division
import numpy as np
from Interface_Plasmon_Simulation.materials import vac, Al
from Interface_Plasmon_Simulation.microscope import Microscope
from Interface_Plasmon_Simulation.bulkMode import bulk_NR
from Interface_Plasmon_Simulation.twoSlabParallel import twoSlabParallel_NR, twoSlabParallel, dispersion
import matplotlib.pyplot as plt
import matplotlib.colors as pltc

#%%
c = 3E8 #[m/s]
hbar=6.582E-16#[eV s]

#%%
x = np.array([-400, 0, 400])*1E-10#np.arange(-5,0,0.5)*1E-10 #[m]
E = np.arange(0.001,17,0.03)
q_perpendicular = np.arange(0,200,0.05)*1E6

#%%
materials = [vac,Al]
for material in materials:
    print(material.name)
    material.set_Ep()
    material.set_eps(E=E)

#%%
microscope = Microscope()
microscope.print_parameters()

#%%
pos_index= 1

#%%
interface, begrenzung = twoSlabParallel_NR(microscope, material=materials,
                   x=x, q_perpendicular=q_perpendicular, E=E)
bulk = bulk_NR(microscope, material=materials,
                   x=x, q_perpendicular=q_perpendicular, E=E)

fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2)
fig.suptitle(r'pos: {} [$\AA$]'.format(x[pos_index]))

img1 = ax1.imshow(interface[:,pos_index,:], aspect='auto', norm=pltc.LogNorm(), origin='lower',
                 extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
ax1.plot(q_perpendicular*c*hbar, q_perpendicular, color='red')
ax1.set_title('Interface')
plt.xlabel('E [eV]')
ax1.set_xlim(0,np.amax(E))
plt.ylabel(r'$q_{y} [m^-]$')
fig.colorbar(img1, ax=ax1)
#    plt.colorbar()

img2 = ax2.imshow(begrenzung[:,pos_index,:], aspect='auto', origin='lower',
             extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
ax2.set_title('Begrenzung')
plt.xlabel('E [eV]')
plt.ylabel(r'$q_{y} [m^-]$')
fig.colorbar(img2, ax=ax2)

img3 = ax3.imshow(bulk[:,pos_index,:], aspect='auto', origin='lower',
             extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
ax3.set_title('Bulk')
plt.xlabel('E [eV]')
plt.ylabel(r'$q_{y} [m^-]$')
fig.colorbar(img3, ax=ax3)

img4 = ax4.imshow(bulk[:,pos_index,:]+begrenzung[:,pos_index,:]+interface[:,pos_index,:], aspect='auto', origin='lower',
             extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
ax4.plot(q_perpendicular*c*hbar, q_perpendicular, color='red')
ax4.set_title('Total')
plt.xlabel('E [eV]')
ax4.set_xlim(0,np.amax(E))
plt.ylabel(r'$q_{y} [m^-]$')
fig.colorbar(img4, ax=ax4)

#%%
interface_s, interface_p = twoSlabParallel(microscope, material=materials,
                   x=x, q_perpendicular=q_perpendicular, E=E)

fig, [ax1, ax2] = plt.subplots(2, 1)
img1 = ax1.imshow(interface_s[:,-1,:], aspect='auto', origin='lower',
                 extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
plt.xlabel('E [eV]')
plt.ylabel(r'$q_{y} [m^-]$')
fig.colorbar(img1, ax=ax1)
#    plt.colorbar()

img2 = ax2.imshow(interface_p[:,-1,:], aspect='auto', origin='lower',
             extent=(np.amin(E),np.amax(E), np.amin(q_perpendicular), np.amax(q_perpendicular)))
plt.xlabel('E [eV]')
plt.ylabel(r'$q_{y} [m^-]$')
fig.colorbar(img2, ax=ax2)
ax2.plot(q_perpendicular*c*hbar, q_perpendicular, color='red')
ax2.set_xlim(0,np.amax(E))

fig2 = plt.figure('dispersion')
plt.plot(E,dispersion(material=materials, E=E))
plt.plot(q_perpendicular*c*hbar, q_perpendicular, color='red')
plt.xlim(0,np.amax(E))
plt.xlabel('E [eV]')
plt.ylabel(r'$q_{||} [m^-]$')

plt.show()