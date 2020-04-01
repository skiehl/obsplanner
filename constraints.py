# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Classes for observational constraints.
"""

from abc import ABCMeta, abstractmethod
from astropy.coordinates import get_sun, get_moon
import astropy.units as u
import numpy as np

import utilities as ut

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

class Constraints(object):
    """List of constraints that are getting jointly evaluated.
    """

    #--------------------------------------------------------------------------
    def __init__(self):
        """Create a Constraints instance.
        """

        self.constraints = []

    #--------------------------------------------------------------------------
    def add(self, constraint):
        """Add a new constraint.
        """

        # TODO: check if Constraint instance (need ABC module)

        self.constraints.append(constraint)

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """Evaluate all constraints jointly.
        
        Parameters
        -----
        source_coord : astropy.SkyCoord
            Coordinates of the target sources.
        telescope : Telescope
            Provides the telescope position and current date and time.
        
        Returns
        -----
        out : numpy.ndarray
            Array of boolean values. One entry for each input coordinate.
            The entry is True, if the source is observable according to all
            constraints, and False otherwise.
        """

        observable = np.logical_and.reduce(
                [constraint.get(source_coord, telescope) \
                 for constraint in self.constraints])

        return observable

#==============================================================================

class Constraint(object, metaclass=ABCMeta):
    """Observational constraint. Defines whether a source is observable from a
    specified location at a specified time or not.
    """

    #--------------------------------------------------------------------------
    @abstractmethod
    def __init__(self):
        """Create Constraint instance.
        
        Notes
        -----
        Abstract method. Parameters depend on the specific constraint.
        """

    #--------------------------------------------------------------------------
    @abstractmethod
    def get(self, source_coord, telescope):
        """Evaluate the constraint for given targets, a specific location and
        time.
        
        Parameters
        -----
        source_coord : astropy.SkyCoord
            Coordinates of the target sources.
        telescope : Telescope
            Provides the telescope position and current date and time.
        
        Returns
        -----
        out : numpy.ndarray
            Array of boolean values. One entry for each input coordinate.
            The entry is True, if the source is observable according to the
            constraint, and False otherwise.
        
        Notes
        -----
        Abstract method. The implementation needs to specified in the
        sub-classes.
        """

#==============================================================================

class ElevationLimit(Constraint):
    """Elevation limit: only sources above a specified elevation are
    observable.
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit):
        """Create ElevationLimit instance.
        
        Parameters
        -----
        limit : float
            Lower elevation limit in degrees.
        """

        self.limit = limit * u.deg

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """Evaluate the constraint for given targets, a specific location and
        time.
        
        Parameters
        -----
        source_coord : astropy.SkyCoord
            Coordinates of the target sources.
        telescope : Telescope
            Provides the telescope position and current date and time.
        
        Returns
        -----
        out : numpy.ndarray
            Array of boolean values. One entry for each input coordinate.
            The entry is True, if the source is observable according to the
            constraint, and False otherwise.
        """

        altaz = source_coord.transform_to(telescope.frame)
        observable = altaz.alt >= self.limit
        #observable = np.array(observable, dtype=float)

        return observable

#==============================================================================

class AirmassLimit(Constraint):
    """Airmass limit: only sources below a specified airmass are observable.
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit, conversion="secz"):
        """Create AirmassLimit instance.
        
        Parameters
        -----
        limit : float
            Upper airmass limit
        conversion : str, default="secz"
            Select the conversion method from altitude to airmass.
            Options: "secz" (default), "Rosenberg", "KastenYoung",
            "Young". See alt_to_airmass() docstring in utilities.py for
            details.
        """

        self.limit = limit
        self.conversion = conversion

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """Evaluate the constraint for given targets, a specific location and
        time.
        
        Parameters
        -----
        source_coord : astropy.SkyCoord
            Coordinates of the target sources.
        telescope : Telescope
            Provides the telescope position and current date and time.
        
        Returns
        -----
        out : numpy.ndarray
            Array of boolean values. One entry for each input coordinate.
            The entry is True, if the source is observable according to the
            constraint, and False otherwise.
        """

        altaz = source_coord.transform_to(telescope.frame)
        observable = ut.alt_to_airmass(altaz.al, conversion=self.conversion) \
                <= self.limit
        #observable = np.array(observable, dtype=float)

        return observable

#==============================================================================

class SunDistance(Constraint):
    """Only sources sufficiently separated from the Sun are observable.
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit):
        """Create SunDistance instance.
        
        Parameters
        -----
        limit : float
            Minimum required angular separation between targets and the Sun in
            degrees.
        """

        self.limit = limit * u.deg

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """Evaluate the constraint for given targets, a specific location and
        time.
        
        Parameters
        -----
        source_coord : astropy.SkyCoord
            Coordinates of the target sources.
        telescope : Telescope
            Provides the telescope position and current date and time.
        
        Returns
        -----
        out : numpy.ndarray
            Array of boolean values. One entry for each input coordinate.
            The entry is True, if the source is observable according to the
            constraint, and False otherwise.
        """

        source_altaz = source_coord.transform_to(telescope.frame)
        sun_altaz = get_sun(telescope.time).transform_to(telescope.frame)
        separation = source_altaz.separation(sun_altaz)
        observable = separation > self.limit
        #observable = np.array(observable, dtype=float)

        return observable

#==============================================================================

class MoonDistance(Constraint):
    """Only sources sufficiently separated from the Moon are observable.
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit):
        """Create MoonDistance instance.
        
        Parameters
        -----
        limit : float
            Minimum required angular separation between targets and the Moon in
            degrees.
        """

        self.limit = limit * u.deg

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """Evaluate the constraint for given targets, a specific location and
        time.
        
        Parameters
        -----
        source_coord : astropy.SkyCoord
            Coordinates of the target sources.
        telescope : Telescope
            Provides the telescope position and current date and time.
        
        Returns
        -----
        out : numpy.ndarray
            Array of boolean values. One entry for each input coordinate.
            The entry is True, if the source is observable according to the
            constraint, and False otherwise.
        """

        source_altaz = source_coord.transform_to(telescope.frame)
        moon_altaz = get_moon(telescope.time).transform_to(telescope.frame)
        separation = source_altaz.separation(moon_altaz)
        observable = separation > self.limit
        #observable = np.array(observable, dtype=float)

        return observable

#==============================================================================

class MoonPolarization(Constraint):
    """Avoid polarized, scattered Moon light. Target sources in a specified
    angular range around 90 degrees separation from the Moon are not
    observable.
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit):
        """Create MoonPolarization instance.
        
        Parameters
        -----
        limit : float
            Angular range to avoid polarized, scattered Moon light. Sources
            within the range (90-limit, 90+limit) degrees separation from the
            Moon are not observable.
        """

        self.limit = limit * u.deg
        self.limit_lo = (90. - limit) * u.deg
        self.limit_hi = (90. + limit) * u.deg

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """Evaluate the constraint for given targets, a specific location and
        time.
        
        Parameters
        -----
        source_coord : astropy.SkyCoord
            Coordinates of the target sources.
        telescope : Telescope
            Provides the telescope position and current date and time.
        
        Returns
        -----
        out : numpy.ndarray
            Array of boolean values. One entry for each input coordinate.
            The entry is True, if the source is observable according to the
            constrains, and False otherwise.
        """

        source_altaz = source_coord.transform_to(telescope.frame)
        moon_altaz = get_moon(telescope.time).transform_to(telescope.frame)
        separation = source_altaz.separation(moon_altaz)
        observable = np.logical_or(
                separation <= self.limit_lo, separation >= self.limit_hi)
        #observable = np.array(observable, dtype=float)

        return observable
