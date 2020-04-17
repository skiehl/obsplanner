# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Classes related to the motion of a telescope.
"""

from abc import ABCMeta, abstractmethod
from astropy.coordinates import AltAz, Angle, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
import astropy.units as u
import numpy as np

# alternative download location since
# http://maia.usno.navy.mil/ser7/finals2000A.all is unavailable:
from astropy.utils import iers
iers.conf.iers_auto_url = \
    'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'

__author__ = "Sebastian Kiehlmann"
__credits__ = ["Sebastian Kiehlmann"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sebastian Kiehlmann"
__email__ = "skiehlmann@mail.de"
__status__ = "Production"

#==============================================================================
# CLASSES
#==============================================================================

class Position(object):
    """Position of a telescope.
    """

    #--------------------------------------------------------------------------
    def __init__(self, coord, loc, time):
        """Create a Position instance.

        Parameters
        -----
        coord : astropy.SkyCoord
            Coordinates of a source the telescope is currently pointed at.
        loc : astropy.EarthLocation
            Location of the telescope on Earth.
        time : astropy.Time
            Current date and time.

        Returns
        -----
        None
        """

        frame = AltAz(obstime=time, location=loc)
        altaz = coord.transform_to(frame)
        self.ra = coord.ra
        self.dec = coord.dec
        self.alt = altaz.alt
        self.az = altaz.az
        self. ha = self.ra - time.sidereal_time('apparent')
        self.za = self.dec - loc.lat

    #--------------------------------------------------------------------------
    def __str__(self):
        """Write out the position in various coordinate systems.

        Returns
        -----
        out : str
        """

        text = 'Right Ascension: {0:s}\n'.format(str(self.ra))
        text += 'Declination:     {0:s}\n'.format(str(self.dec))
        text += 'Altitude:       {0:6.2f}\n'.format(self.alt)
        text += 'Azimuth:        {0:6.2f}\n'.format(self.az)
        text += 'Hour angle:     {0:6.2f} h\n'.format(self.ha.hour)
        text += 'Zenith angle:   {0:6.2f}'.format(self.za)

        return text

#==============================================================================

class Telescope(object, metaclass=ABCMeta):
    """Telescope base class.

    Provides general and abstract methods for Telescope instances.
    Child classes for different telescope mounts need to specify the abstract
    methods.
    """

    #--------------------------------------------------------------------------
    @abstractmethod
    def __init__(
            self, lat, lon, height, name=''):
        """Create a Telescope instance.

        Parameters
        -----
        lat : str or astropy.Angle
            Latitude of telescope location. String input needs to be consistent
            with astropy.Angle definition.
        lon : str or astropy.Angle
            Longitude of telescope location. String input needs to be
            consistent with astropy.Angle definition.
        height : float
            Height of telescope location in meters.
        name : str, default=''
            Name of the telescope/observatory.

        Returns
        -----
        None

        Notes
        -----
        Abstract method. Mount specific attributes are defined in the child
        classes.
        """

        lat = Angle(lat)
        lon = Angle(lon)
        height = height * u.m

        self.loc = EarthLocation(
                lat=lat, lon=lon, height=height)
        self.name = name
        self.dome = False
        self.time = None
        self.frame = None
        self.pos = None

        self.slew_model_dome = None

        print('Telescope: {0:s} created.'.format(self.name))

    #--------------------------------------------------------------------------
    def __str__(self):
        """Write out location and current position of the telescope.

        Returns
        -----
        out : str
        """

        text = 'Telescope:  {0:s}\n'.format(self.name)
        text += 'Longitude:  {0:s}\n'.format(
                self.loc.lon.to_string(decimal=False))
        text += 'Longitude:  {0:s}\n'.format(
                self.loc.lat.to_string(decimal=False))
        text += 'Height:     {0:.2f}\n'.format(self.loc.height)
        text += 'Mount:      {0:s}\n'.format(self.mount)

        if self.time is not None:
            text += '\nTime:       {0:s}\n'.format(str(self.time))

        if self.pos is not None:
            text += '\n' + self.pos.__str__()

        return text

    #--------------------------------------------------------------------------
    @abstractmethod
    def set_slew_model(self, axis, slew_model):
        """Set a slew model for a specific axis.
        This slew model will define how to calculate the time needed to slew
        from one position to another.

        Parameters
        -----
        axis : str
            Specify axis: 'ra', 'dec', 'dome'
        slew_model : SlewModel
            Slew model used for the specified axis.

        Returns
        -----
        None

        Notes
        -----
        Abstract method. Axis is mount specific.
        """

        if axis == 'dome':
            self.dome = True
            self.slew_model_dome = slew_model
            print('Telescope: Dome slew model set.')

        elif axis == 'ra':
            # implement in child class
            pass

        elif axis == 'dec':
            # implement in child class
            pass

        else:
            raise ValueError(
                    "Unknown axis: '{0:s}'".format(axis))

        return None

    #--------------------------------------------------------------------------
    @abstractmethod
    def is_set(self):
        """Check if slew models have been set.

        Returns
        -----
        out : bool
            True, if slew model has been set for each axis. False, otherwise.

        Notes
        -----
        Abstract method. Mount dependent specifics need to be defined in child
        classes.
        """

        return False

    #--------------------------------------------------------------------------
    def set_time(self, time):
        """Set current time.

        Parameters
        -----
        time : astropy.Time
            Date and time in UTC. Accepts any input format that astropy.Time is
            accepting.

        Returns
        -----
        None
        """

        self.time = Time(time, scale='utc', location=self.loc)
        self.frame = AltAz(obstime=self.time, location=self.loc)

    #--------------------------------------------------------------------------
    def set_pos(self, coord):
        """Set current position of telescope.

        Parameters
        -----
        coord : astropy.SkyCoord
            Coordinates of a source the telescope is currently pointed at.

        Returns
        -----
        None
        """

        if self.time is None:
            raise ValueError(
                    "Set time first before setting a position.")

        self.pos = Position(coord, self.loc, self.time)

        return None

    #--------------------------------------------------------------------------
    @abstractmethod
    def set_to_zenith(self):
        """Set telescope position to zenith.

        Returns
        -----
        None

        Notes
        -----
        Abstract method. Mount dependent specifics need to be defined in child
        classes.
        """

        # implement in child class
        pass

    #--------------------------------------------------------------------------
    @abstractmethod
    def get_slew_time(self, coord):
        """Calculate slew time from current position to given coordinates.

        Parameters
        -----
        coord : astropy.SkyCoord
            Source coordinates.

        Returns
        -----
        out : numpy.ndarray
            Slew time to each source from the current position.

        Notes
        -----
        Abstract method. Mount dependent specifics need to be defined in child
        classes.
        """

        # implement in child class
        pass

    #--------------------------------------------------------------------------
    def get_time(self):
        """Get the telescope's current date and time.

        Returns
        -----
        out : astropy.Time
            Current date and time of the telescope.
        """

        if self.time is None:
            raise ValueError('No date and time set yet.')

        return self.time

    #--------------------------------------------------------------------------
    def next_sun_set_rise(self, twilight):
        """
        """

        if isinstance(twilight, float):
            sunset = twilight * u.deg
        elif twilight == 'astronomical':
            sunset = -18. * u.deg
        elif twilight == 'nautical':
            sunset = -12. * u.deg
        elif twilight == 'civil':
            sunset = -6. * u.deg
        elif twilight == 'sunset':
            sunset = 0. * u.deg
        else:
            raise ValueError(
                "Either set a float or chose from 'astronomical', 'nautical'" \
                " 'civil', or 'sunset'.")

        print('Current time:    ', self.time)
        time = self.time + np.arange(0., 48., 0.1) * u.h
        frame = AltAz(obstime=time, location=self.loc)
        sun = get_sun(time).transform_to(frame)
        sun_alt = sun.alt
        night = sun_alt < sunset

        if np.all(night):
            print('NOTE: Sun never sets!')
            return False

        if np.all(~night):
            print('WARNING: Sun never rises!')
            return False

        # current time is at night, find next sun set and following sun rise:
        if night[0]:
            print('NOTE: current time is night time.')
            # index of next sun rise:
            i = np.argmax(~night)
            # index of next sun set:
            i += np.argmax(night[i:])
            # index of following sun rise:
            j = i + np.argmax(~night[i:]) - 1

        # current time is at day, find next sun set and sun rise:
        else:
            # index of next sun set:
            i = np.argmax(night)
            # index of following sun rise:
            j = i + np.argmax(~night[i:]) - 1


        # interpolate linearly:
        interp = (sunset.value - sun_alt[i-1].deg) \
                / (sun_alt[i].deg - sun_alt[i-1].deg)
        time_sunset = time[i-1] + (time[i] - time[i-1]) * interp
        interp = (sunset.value - sun_alt[j].deg) \
                / (sun_alt[j+1].deg - sun_alt[j].deg)
        time_sunrise = time[j] + (time[j+1] - time[j]) * interp

        print('Next night start:', time_sunset)
        print('Next night stop: ', time_sunrise)

        return time_sunset, time_sunrise

#==============================================================================

class TelescopeEq(Telescope):
    """Telescope with equatorial mount.
    """

    mount = 'equatorial'

    #--------------------------------------------------------------------------
    def __init__(
            self, lat, lon, height, name=''):
        """Create a Telescope instance.

        Parameters
        -----
        lat : str or astropy.Angle
            Latitude of telescope location. String input needs to be consistent
            with astropy.Angle definition.
        lon : str or astropy.Angle
            Longitude of telescope location. String input needs to be
            consistent with astropy.Angle definition.
        height : float
            Height of telescope location in meters.
        name : str, default=''
            Name of the telescope/observatory.

        Returns
        -----
        None
        """

        super().__init__(lat, lon, height, name=name)
        self.slew_model_ra = None
        self.slew_model_dec = None

    #--------------------------------------------------------------------------
    def set_slew_model(self, axis, slew_model):
        """Set a slew model for a specific axis.
        This slew model will define how to calculate the time needed to slew
        from one position to another.

        Parameters
        -----
        axis : str
            Specify axis: 'ra', 'dec', 'dome'
        slew_model : SlewModel
            Slew model used for the specified axis.

        Returns
        -----
        None
        """

        if axis == 'dome':
            self.dome = True
            self.slew_model_dome = slew_model
            print('Telescope: Dome slew model set.')

        elif axis == 'ra':
            self.slew_model_ra = slew_model
            print('Telescope: RA slew model set.')

        elif axis == 'dec':
            self.slew_model_dec = slew_model
            print('Telescope: Dec slew model set.')

        else:
            raise ValueError(
                    "Unsupported axis '{0:s}' for mount '{1:s}'.". format(
                            axis, self.mount))

        return None

    #--------------------------------------------------------------------------
    def is_set(self):
        """Check if slew models have been set.

        Returns
        -----
        out : bool
            True, if slew model has been set for each axis. False, otherwise.
        """

        ready = True

        if self.slew_model_ra is None:
            print('WARNING: No slew model for RA axis set.')
            ready = False

        if self.slew_model_dec is None:
            print('WARNING: No slew model for DEC axis set.')
            ready = False

        if self.slew_model_dome is None:
            print('WARNING: No slew model for dome set.')
            # only print(warning, dome is optional

        return ready

    #--------------------------------------------------------------------------
    def set_to_zenith(self):
        """Set telescope position to zenith.

        Returns
        -----
        None
        """

        if self.time is None:
            raise ValueError(
                    "Set time first before setting to zenith.")

        ra = self.time.sidereal_time('apparent')
        dec = self.loc.lat
        coord = SkyCoord(ra=ra, dec=dec)
        self.set_pos(coord)

    #--------------------------------------------------------------------------
    def get_slew_time(self, coord):
        """Calculate slew time from current position to given coordinates.

        Parameters
        -----
        coord : astropy.SkyCoord
            Source coordinates.

        Returns
        -----
        out : numpy.ndarray
            Slew time to each source from the current position.
        """

        # calculate slew times of individual axes:
        rot_ra = coord.ra - self.pos.ra
        rot_dec = self.pos.dec - coord.dec
        t_ra = self.slew_model_ra.get_slew_time(rot_ra)
        t_dec = self.slew_model_dec.get_slew_time(rot_dec)
        time = [t_ra, t_dec]

        if self.dome:
            t_dome = self.slew_model_dome.get_slew_time(rot_ra)
            time.append(t_dome)

        # take maximum of all axes:
        time = np.maximum.reduce(time)

        return time

#==============================================================================

class SlewModel(object, metaclass=ABCMeta):
    """Slew model base class.
    """

    #--------------------------------------------------------------------------
    @abstractmethod
    def __init__(self):
        """Create a SlewModel instance.
        """

        return None

    #--------------------------------------------------------------------------
    @abstractmethod
    def __str__(self):
        """Write out information about the slew model.

        Returns
        -----
        out : str
        """

        return 'SlewModel instance.'

    #--------------------------------------------------------------------------
    @abstractmethod
    def get_slew_time(self, angle):
        """Calculate the slew time.

        Parameters
        -----
        angle : astropy.Angle
            Angle(s) by which to slew.

        Returns
        -----
        out : numpy.ndarray
            Time needed to slew by the given angle(s).
        """

        return None


#==============================================================================

class SlewModelLinear(SlewModel):
    """Linear slew model.
    """

    #--------------------------------------------------------------------------
    def __init__(self, init_p, vel_p, init_n=None, vel_n=None):
        """Create SlewModelLinear instance.

        Parameters
        -----
        init_p : float
            Initial time in seconds needed to start slewing in positive
            direction.
        vel_p : float
            Slew velocity in degrees per second for slewing in positive
            direction.
        init_n : float, default=None
            Initial time in seconds needed to start slewing in negaitve
            direction. Same as positive direction, if not set.
        vel_n : float, default=None
            Slew velocity in degrees per second for slewing in negative
            direction. Same as positive direction, if not set.
        """

        self.init_p = init_p * u.s
        self.vel_p = vel_p * u.deg / u.s
        if init_n is None:
            self.init_n = init_p  * u.s
        else:
            self.init_n = init_n  * u.s
        if vel_n is None:
            self.vel_n = vel_p * u.deg / u.s
        else:
            self.vel_n = vel_n * u.deg / u.s

        return None

    #--------------------------------------------------------------------------
    def __str__(self):
        """Write out information about the slew model.

        Returns
        -----
        out : str
        """

        text = 'SlewModel: linear\n'
        if self.vel_n == self.vel_p and self.init_n == self.init_p:
            text += 'Velocity:        {0:10.4f}\n'.format(
                    self.vel_p)
            text += 'Additional time: {0:10.4f}\n'.format(self.init_p)
        else:
            text += 'Positive direction:\n'
            text += 'Velocity:        {0:10.4f}\n'.format(
                    self.vel_p)
            text += 'Additional time: {0:10.4f}\n'.format(self.init_p)
            text += 'Negative direction:\n'
            text += 'Velocity:        {0:10.4f}\n'.format(
                    self.vel_n)
            text += 'Additional time: {0:10.4f}\n'.format(self.init_n)

        return text

    #--------------------------------------------------------------------------
    def get_slew_time(self, angle):
        """Calculate the slew time.

        Linear model:
        .. math::
            t = t_0 + a / v,
        where t_0 is an initial time offset, a is the angle to rotate by, and
        v is the rotation velocity. Parameters can be set individually for
        rotations in positive and negative directions.

        Parameters
        -----
        angle : astropy.Angle
            Angle(s) by which to slew.

        Returns
        -----
        out : numpy.ndarray
            Time needed to slew by the given angle(s).
        """

        time = np.ones(angle.size) * u.s

        if angle.size == 1:
            if angle > 0:
                time = angle  / self.vel_p + self.init_p
            elif angle < 0:
                time = angle  / self.vel_n + self.init_n
            else:
                time = 0 * u.s

        elif angle.size > 1:
            sel = angle == 0
            time[sel] = 0

            sel = angle > 0
            time[sel] = angle[sel] / self.vel_p + self.init_p

            sel = angle < 0
            time[sel] = -angle[sel] / self.vel_n + self.init_n

        return time