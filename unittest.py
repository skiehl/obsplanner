# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Module description.
"""

import numpy as np

from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

import constraints as c
import sources as s
import telescope as t
import utilities as ut
import astropy.units as u


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

#------------------------------------------------------------------------------
print('--------------------------\nTest: utilities.py\n')

za = 60. *u.deg
print('Zenith angle {0:.1f} degrees corresponds to airmass:'.format(za.value))
for conv in ['secz', 'Rosenberg', 'KastenYoung', 'Young']:
    airmass = ut.za_to_airmass(za, conversion=conv)
    print('{0:6.2f}  ({1:s})'.format(airmass, conv))

#------------------------------------------------------------------------------
print('--------------------------\nTest: sources.py\n')

dtype = [('name', 'S30'), ('ra', float), ('dec', float), ('expt', float),
         ('expn', int)]
targets = np.loadtxt(
        'sourcelists/targets.csv', dtype=dtype, skiprows=1, delimiter=',',
        usecols=range(5))

sources = s.Sources(targets['name'], targets['ra'], targets['dec'])

#------------------------------------------------------------------------------
print('--------------------------\nTest: telescope.py\n')

print('Testing class: Position\n')
coord = SkyCoord('75d12m14.1035s', '24d53m57s')
loc = EarthLocation(lat='24:53:57 deg', lon='35:12:43 deg', height=1750)
time = Time('2019-11-01 00:00:00', location=loc)
pos = t.Position(coord, loc, time)
print(pos)

print('\nTesting class: Telescope\n')
telescope = t.TelescopeEq(
        '24:53:57 deg', '35:12:43 deg', 1750, name='Skinakas')
telescope.set_time('2019-11-01 00:00:00')
telescope.set_to_zenith()
print(telescope)

print('Testing class: SlewModelLinear')
slewmodel = t.SlewModelLinear(18., 0.5)
print(slewmodel)
slewmodel = t.SlewModelLinear(18., 0.5, 18.1, 0.6)
print(slewmodel)

print('Angle        Slew time')
angle = sources.coord.ra[0] - sources.coord.ra[1]
print(angle, slewmodel.get_slew_time(angle))
angle = sources.coord.ra[3] - sources.coord.ra[0:5]
for a, st in zip(angle, slewmodel.get_slew_time(angle)):
    print(a, st)

print('\nTest: adding slew model to telescope')
slewmodel_ra = t.SlewModelLinear(2., 0.5, 2.1, 0.6)
slewmodel_dec = t.SlewModelLinear(1.9, 0.4, 1.8, 0.5)
slewmodel_dome = t.SlewModelLinear(20., 0.1)
telescope.set_slew_model('ra', slewmodel_ra)
telescope.set_slew_model('dec', slewmodel_dec)
#telescope.set_slew_model('x', slewmodel_dec)
telescope.set_slew_model('dome', slewmodel_dome)
print(telescope.is_set())


print(telescope.get_slew_time(sources.coord))

#------------------------------------------------------------------------------
print('--------------------------\nTest: constraints.py\n')

elevation_limit = c.ElevationLimit(30.)
print(elevation_limit.get(sources.coord, telescope))
airmass_limit = c.AirmassLimit(2.)
print(airmass_limit.get(sources.coord, telescope))
sun_distance = c.SunDistance(20.)
print(sun_distance.get(sources.coord, telescope))
moon_distance = c.MoonDistance(20.)
print(moon_distance.get(sources.coord, telescope))
moon_polarization = c.MoonPolarization(10.)
print(moon_polarization.get(sources.coord, telescope))

constraints = c.Constraints()
constraints.add(airmass_limit)
constraints.add(sun_distance)
constraints.add(moon_distance)
constraints.add(moon_polarization)
print(constraints.get(sources.coord, telescope))

