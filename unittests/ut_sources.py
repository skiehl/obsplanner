# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Module description.
"""

import sys
sys.path.append('../')
import astropy.units as u
import numpy as np

import sources as s

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

print('\nTEST sources.Sources:')

dtype = [('name', 'S30'), ('ra', float), ('dec', float), ('expt', float),
         ('expn', int)]
targets = np.loadtxt(
        'sourcelists/targets.csv', dtype=dtype, skiprows=1, delimiter=',',
        usecols=range(5))

print('Initialize instance:')
sources = s.Sources(
    targets['name'], targets['ra'], targets['dec'], targets['expt'],
    targets['expn'])

print('Length i.e. size of instance:')
print(len(sources))
print(sources.size)
