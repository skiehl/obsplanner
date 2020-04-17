# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Module description.
"""

import sys
sys.path.append('../')
import astropy.units as u

import utilities as ut

__author__ = "Sebastian Kiehlmann"
__credits__ = ["Sebastian Kiehlmann"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sebastian Kiehlmann"
__email__ = "skiehlmann@mail.de"
__status__ = "Production"

#==============================================================================
# MAIN
#==============================================================================

print('\nTEST utilities.za_to_airmass():')

za = 60. * u.deg
print('Zenith angle {0:.1f} degrees corresponds to airmass:'.format(za.value))

for conv in ['secz', 'Rosenberg', 'KastenYoung', 'Young']:
    airmass = ut.za_to_airmass(za, conversion=conv)
    print('{0:6.2f}  ({1:s})'.format(airmass, conv))

#------------------------------------------------------------------------------

print('\nTEST utilities.alt_to_airmass():')

za = 30. * u.deg
print('Altitude {0:.1f} degrees corresponds to airmass:'.format(za.value))

for conv in ['secz', 'Rosenberg', 'KastenYoung', 'Young']:
    airmass = ut.alt_to_airmass(za, conversion=conv)
    print('{0:6.2f}  ({1:s})'.format(airmass, conv))