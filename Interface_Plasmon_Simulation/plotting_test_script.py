# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 18:01:15 2017

@author: Eric Hoglund
"""

import numpy as np
from materials import Al, GB, vac
from microscope import microscope
from Bulk_plasmon import bulk_plasmon_double_differential_cross_section
from normal_incidence_wRetardation import double_differential_cross_section_normalIncidence


scope = microscope(keV=100)
scope.print_parameters()

class create_spectrum():
    def __init__(self, q_perp_max=2E10, E_min=5, E_max=20):
        self.q_perp = np.linspace(-1*q_perp_max, q_perp_max,800)[:,None]#microscope.k0*microscope.collection_angle
        self.E = np.linspace(E_min,E_max,400)

dimensions = create_spectrum(q_perp_max=2E10, E_min=5, E_max=20)
print('theta_max: %.4f'%(dimensions.q_perp.max()/scope.k0))
print('(111) @',2*np.pi/2.338,'[1/A]')
print()

#MR_Spectra = bulk_plasmon_double_differential_cross_section(microscope(),create_spectrum(),Al)
MR_Spectra = double_differential_cross_section_normalIncidence(scope, dimensions, [vac,Al], 20*10E-9)

import matplotlib.pyplot as plt
def plot_limit_lines():
    c=3E8#[m/s]
    e=1.602E-19 #[C]
    m0=9.11E-31#[kg]
    plt.plot(np.linspace(dimensions.E.min(), dimensions.E.max(),dimensions.E.size()), np.linspace(dimensions.E.min(), dimensions.E.max(),dimensions.E.size())/(m0*scope.v*c/e)*scope.k0 *10**-10, color='b')


fig, ax = plt.subplots()
#print(MR_Spectra.eps.real)
plt.plot(MR_Spectra.eps[1][0].real, label=r'$\epsilon_r q=0$')
plt.plot(MR_Spectra.eps[1][200].real, label=r'$\epsilon_r q=200$')
plt.plot(MR_Spectra.eps[1][0].conj().imag, label=r'$\epsilon_i q=0$')
plt.plot(MR_Spectra.eps[1][200].conj().imag, label=r'$\epsilon_i q=200$')
plt.legend()


bounds = [dimensions.E.min(), dimensions.E.max(), dimensions.q_perp.min() *10**-10, dimensions.q_perp.max() *10**-10]

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



######################################
######################################
############   plot Ep(q)  ###########
######################################
######################################
# check that it goes up to ~20-30  ###

plt.show()