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


scope = microscope()

class create_spectrum():
    def __init__(self):
        self.q_perp = np.linspace(-3E10,3E10,400)[:,None]#microscope.k0*microscope.collection_angle
        self.E = np.linspace(10,20,100)

#MR_Spectra = bulk_plasmon_double_differential_cross_section(microscope(),create_spectrum(),Al)
MR_Spectra = double_differential_cross_section_normalIncidence(microscope(),create_spectrum(),[vac,Al],300E9)
print(vac.E_fermi)
import matplotlib.pyplot as plt
fig, ax= plt.subplots()
plt.imshow(MR_Spectra.surface_mode(),origin='lower', aspect='auto', cmap=plt.get_cmap('hot'),
           extent=[10,20,0,4E10*10**-10])
plt.colorbar(label=r'$\frac{dP^3}{dz dE dq_y} \left[ \frac{1}{eV nm}\right] $')
plt.title('Bulk')
plt.ylabel(r'$q_y [A^-]$')
plt.xlabel(r'$E [eV]$')
plt.show()