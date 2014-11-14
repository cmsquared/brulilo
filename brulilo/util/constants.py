"""
Some simple constants that we need.  Should probably use something more formal
here...
"""
import numpy as np

MeV2eV = 1e6
eV2erg = 1.6021766e-12
MeV2erg = MeV2eV * eV2erg

# cgs units
electron_mass = 9.10938215e-28     # g
amu = 1.660538921e-24              # g

light_speed = 2.99792458e10       # cm / s
planck = 6.62606957e-27            # erg * s
planck_bar = planck / (2*np.pi)    # erg * s
boltzmann = 1.3806488e-16          # erg / K

avogadro = 6.0221413e23
