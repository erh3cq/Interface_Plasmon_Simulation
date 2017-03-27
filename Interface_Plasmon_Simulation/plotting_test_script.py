# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 18:01:15 2017

@author: Eric Hoglund
"""

import numpy as np
from materials import Al, GB, vac
from microscope import microscope
from Bulk_Plasmon import bulk_plasmon_double_differential_cross_section
from normal_incidence_wRetardation import double_differential_cross_section_normalIncidence


scope = microscope()

q_bound = 4E10
class create_spectrum():
    def __init__(self):
        self.q_perp = np.linspace(-1*q_bound, q_bound,400)[:,None]#microscope.k0*microscope.collection_angle
        self.E = np.linspace(1,20,200)

#MR_Spectra = bulk_plasmon_double_differential_cross_section(microscope(),create_spectrum(),Al)
MR_Spectra = double_differential_cross_section_normalIncidence(microscope(),create_spectrum(),[vac,Al],20*10E-9)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
#print(MR_Spectra.eps.real)
plt.plot(MR_Spectra.eps[1][0].real, label=r'$\epsilon_r q=0$')
plt.plot(MR_Spectra.eps[1][200].real, label=r'$\epsilon_r q=200$')
plt.plot(MR_Spectra.eps[1][0].conj().imag, label=r'$\epsilon_i q=0$')
plt.plot(MR_Spectra.eps[1][200].conj().imag, label=r'$\epsilon_i q=200$')
plt.legend()


bounds = [1,20,-1*q_bound *10**-10,q_bound *10**-10]

fig, ax= plt.subplots(2,2)

im = ax[0,0].imshow(MR_Spectra.bulk_mode(), origin='lower', aspect='auto', cmap=plt.get_cmap('hot'),
           extent=bounds)
fig.colorbar(im, ax=ax[0,0], label=r'$\frac{dP^3}{dz dE dq_y} \left[ eV\ nm^{-1} \right] $')
ax[0,0].set_title('Bulk')
ax[0,0].set_ylabel(r'$q_y [A^-]$')
ax[0,0].set_xlabel(r'$E [eV]$')

im = ax[0,1].imshow(MR_Spectra.surface_mode(), origin='lower', aspect='auto', cmap=plt.get_cmap('hot'),
           extent=bounds)
fig.colorbar(im, ax=ax[0,1], label=r'$\frac{dP^3}{dz dE dq_y} \left[ eV\ nm^{-1} \right] $')
ax[0,1].set_title('Surface')
ax[0,1].set_ylabel(r'$q_y [A^-]$')
ax[0,1].set_xlabel(r'$E [eV]$')

im = ax[1,0].imshow(MR_Spectra.guidedLight1(), origin='lower', aspect='auto', cmap=plt.get_cmap('hot'),
           extent=bounds)
fig.colorbar(im, ax=ax[1,0], label=r'$\frac{dP^3}{dz dE dq_y} \left[ eV\ nm^{-1} \right] $')
ax[1,0].set_title('Guided light 1')
ax[1,0].set_ylabel(r'$q_y [A^-]$')
ax[1,0].set_xlabel(r'$E [eV]$')

im = ax[1,1].imshow(MR_Spectra.guidedLight2(), origin='lower', aspect='auto', cmap=plt.get_cmap('hot'),
           extent=bounds)
fig.colorbar(im, ax=ax[1,1], label=r'$\frac{dP^3}{dz dE dq_y} \left[ eV\ nm^{-1} \right] $')
ax[1,1].set_title('Guided light 1')
ax[1,1].set_ylabel(r'$q_y [A^-]$')
ax[1,1].set_xlabel(r'$E [eV]$')

fig, ax= plt.subplots()
im = ax.imshow(MR_Spectra.total(), origin='lower', aspect='auto', cmap=plt.get_cmap('hot'),
           extent=bounds)
fig.colorbar(im, ax=ax, label=r'$\frac{dP^3}{dz dE dq_y} \left[ eV\ nm^{-1} \right] $')
ax.set_title('Total')
ax.set_ylabel(r'$q_y [A^-]$')
ax.set_xlabel(r'$E [eV]$')

fig, ax= plt.subplots()
im = ax.imshow(np.log10(MR_Spectra.total()), origin='lower', aspect='auto', cmap=plt.get_cmap('hot'),
           extent=bounds)
fig.colorbar(im, ax=ax, label=r'$\frac{dP^3}{dz dE dq_y} \left[ eV\ nm^{-1} \right] $')
ax.set_title('log10 Total')
ax.set_ylabel(r'$q_y [A^-]$')
ax.set_xlabel(r'$E [eV]$')

plt.show()