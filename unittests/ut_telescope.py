# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Module description.
"""

import sys
sys.path.append('../')
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun
from astropy.time import Time
import numpy as np

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

print('\nTEST telescope.Position:')
coord = SkyCoord('75d12m14.1035s', '24d53m57s')
loc = EarthLocation(lat='24:53:57 deg', lon='35:12:43 deg', height=1750)
time = Time('2019-11-01 00:00:00', location=loc)
pos = t.Position(coord, loc, time)
print(pos)

#------------------------------------------------------------------------------

print('\nCreate source list for testing slew models and telescope:')

dtype = [('name', 'S30'), ('ra', float), ('dec', float), ('expt', float),
         ('expn', int)]
targets = np.loadtxt(
        'sourcelists/targets.csv', dtype=dtype, skiprows=1, delimiter=',',
        usecols=range(5))
sources = s.Sources(
    targets['name'], targets['ra'], targets['dec'], targets['expt'],
    targets['expn'])

#------------------------------------------------------------------------------

print('\nTEST telescope.SlewModelLinear:')
print('Initialize slew model with two parameters:')
slewmodel = t.SlewModelLinear(18., 0.5)
print(slewmodel)

print('Initialize slew model with positive and negative parameters:')
slewmodel = t.SlewModelLinear(18., 0.5, 18.1, 0.6)
print(slewmodel)

print('\nTEST telescope.SlewModelLinear.get_slew_time():')
print('Angle        Slew time')
angle = sources.coord.ra[0] - sources.coord.ra[1]
print(angle, slewmodel.get_slew_time(angle))
angle = sources.coord.ra[3] - sources.coord.ra[0:5]
for a, st in zip(angle, slewmodel.get_slew_time(angle)):
    print(a, st)

#------------------------------------------------------------------------------

print('\nTEST telescope.TelescopeEq:')
print('Initialize Telescope instance:')
telescope = t.TelescopeEq(
        '24:53:57 deg', '35:12:43 deg', 1750, name='Skinakas')

print('\nTEST telescope.Telescope.set_time():')
telescope.set_time('2019-11-01 00:00:00')

print('\nTEST telescope.Telescope.set_to_zenith():')
telescope.set_to_zenith()

print('\nTEST Print telescope:')
print(telescope)

print('\nTEST telescope.Telescope.is_set():')
print(telescope.is_set())

print('\nTEST telescope.TelescopeEq.set_slew_model():')
slewmodel_ra = t.SlewModelLinear(2., 0.5, 2.1, 0.6)
slewmodel_dec = t.SlewModelLinear(1.9, 0.4, 1.8, 0.5)
slewmodel_dome = t.SlewModelLinear(20., 0.1)
telescope.set_slew_model('ra', slewmodel_ra)
telescope.set_slew_model('dec', slewmodel_dec)
print(telescope.is_set())

#telescope.set_slew_model('x', slewmodel_dec) # test invalid input

telescope.set_slew_model('dome', slewmodel_dome)
print(telescope.is_set())

print('\nTEST telescope.Telescope.next_sun_set_rise():')
print('Test when current time is during night time:')
twilight = 'astronomical'
print('Twilight:        ', twilight)
print('Current time:    ', telescope.get_time())
sun_set, sun_rise = telescope.next_sun_set_rise(twilight)
frame = AltAz(obstime=sun_set, location=loc)
sun = get_sun(sun_set).transform_to(frame)
print('Time of Sun set: ', sun_set)
print('Sun altitude:    ', sun.alt.deg)
frame = AltAz(obstime=sun_rise, location=loc)
sun = get_sun(sun_rise).transform_to(frame)
print('Time of Sun rise:', sun_set)
print('Sun altitude:    ', sun.alt.deg)

print('\nTest when current time is during day time:')
telescope.set_time('2019-11-01 12:00:00')
print('Twilight:        ', twilight)
print('Current time:    ', telescope.get_time())
sun_set, sun_rise = telescope.next_sun_set_rise(twilight)
frame = AltAz(obstime=sun_set, location=loc)
sun = get_sun(sun_set).transform_to(frame)
print('Time of Sun set: ', sun_set)
print('Sun altitude:    ', sun.alt.deg)
frame = AltAz(obstime=sun_rise, location=loc)
sun = get_sun(sun_rise).transform_to(frame)
print('Time of Sun rise:', sun_set)
print('Sun altitude:    ', sun.alt.deg)

print('\nTest different twilights:')
twilight = 'nautical'
print('Twilight:        ', twilight)
print('Current time:    ', telescope.get_time())
sun_set, sun_rise = telescope.next_sun_set_rise(twilight)
frame = AltAz(obstime=sun_set, location=loc)
sun = get_sun(sun_set).transform_to(frame)
print('Time of Sun set: ', sun_set)
print('Sun altitude:    ', sun.alt.deg)
frame = AltAz(obstime=sun_rise, location=loc)
sun = get_sun(sun_rise).transform_to(frame)
print('Time of Sun rise:', sun_set)
print('Sun altitude:    ', sun.alt.deg)

print('\nTest different twilights:')
twilight = 'civil'
print('Twilight:        ', twilight)
print('Current time:    ', telescope.get_time())
sun_set, sun_rise = telescope.next_sun_set_rise(twilight)
frame = AltAz(obstime=sun_set, location=loc)
sun = get_sun(sun_set).transform_to(frame)
print('Time of Sun set: ', sun_set)
print('Sun altitude:    ', sun.alt.deg)
frame = AltAz(obstime=sun_rise, location=loc)
sun = get_sun(sun_rise).transform_to(frame)
print('Time of Sun rise:', sun_set)
print('Sun altitude:    ', sun.alt.deg)

print('\nTest different twilights:')
twilight = 'sunset'
print('Twilight:        ', twilight)
print('Current time:    ', telescope.get_time())
sun_set, sun_rise = telescope.next_sun_set_rise(twilight)
frame = AltAz(obstime=sun_set, location=loc)
sun = get_sun(sun_set).transform_to(frame)
print('Time of Sun set: ', sun_set)
print('Sun altitude:    ', sun.alt.deg)
frame = AltAz(obstime=sun_rise, location=loc)
sun = get_sun(sun_rise).transform_to(frame)
print('Time of Sun rise:', sun_set)
print('Sun altitude:    ', sun.alt.deg)

print('\nTest different twilights:')
twilight = -4.
print('Twilight:        ', twilight)
print('Current time:    ', telescope.get_time())
sun_set, sun_rise = telescope.next_sun_set_rise(twilight)
frame = AltAz(obstime=sun_set, location=loc)
sun = get_sun(sun_set).transform_to(frame)
print('Time of Sun set: ', sun_set)
print('Sun altitude:    ', sun.alt.deg)
frame = AltAz(obstime=sun_rise, location=loc)
sun = get_sun(sun_rise).transform_to(frame)
print('Time of Sun rise:', sun_set)
print('Sun altitude:    ', sun.alt.deg)

print('\nTEST telescope.TelescopeEq.get_slew_time():')
print(telescope.get_slew_time(sources.coord))

print('\nTEST telescope.Telescope.set_pos():')
print('Set to position of first source')
telescope.set_pos(sources.coord[0])
print('Get slew times again:')
print(telescope.get_slew_time(sources.coord))