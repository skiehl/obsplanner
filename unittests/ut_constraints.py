# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Module description.
"""

import sys
sys.path.append('../')
import numpy as np

import constraints as c
import sources as s
import telescope as t

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

print('\nCreate source list for testing constraints:')

dtype = [('name', 'S30'), ('ra', float), ('dec', float), ('expt', float),
         ('expn', int)]
targets = np.loadtxt(
        'sourcelists/targets.csv', dtype=dtype, skiprows=1, delimiter=',',
        usecols=range(5))
sources = s.Sources(
    targets['name'], targets['ra'], targets['dec'], targets['expt'],
    targets['expn'])

print('Create telescope for testing constraints:')
telescope = t.TelescopeEq(
        '24:53:57 deg', '35:12:43 deg', 1750, name='Skinakas')
telescope.set_time('2019-11-01 00:00:00')

#------------------------------------------------------------------------------

print('\nTEST constraints.ElevationLimit:')
elevation_limit = c.ElevationLimit(30.)
print(elevation_limit.get(sources.coord, telescope))

#------------------------------------------------------------------------------

print('\nTEST constraints.AirmassLimit:')
airmass_limit = c.AirmassLimit(2.)
print(airmass_limit.get(sources.coord, telescope))

#------------------------------------------------------------------------------

print('\nTEST constraints.SunDistance:')
sun_distance = c.SunDistance(20.)
print(sun_distance.get(sources.coord, telescope))

#------------------------------------------------------------------------------

print('\nTEST constraints.MoonDistance:')
moon_distance = c.MoonDistance(20.)
print(moon_distance.get(sources.coord, telescope))

#------------------------------------------------------------------------------

print('\nTEST constraints.MoonPolarization:')
moon_polarization = c.MoonPolarization(10.)
print(moon_polarization.get(sources.coord, telescope))

#------------------------------------------------------------------------------

print('\nTEST constraints.Constraints:')
constraints = c.Constraints()
constraints.add(airmass_limit)
constraints.add(sun_distance)
constraints.add(moon_distance)
constraints.add(moon_polarization)
print(constraints.get(sources.coord, telescope))
