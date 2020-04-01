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
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self):
        """
        """

        self.constraints = []

    #--------------------------------------------------------------------------
    def add(self, constraint):
        """
        """

        # TODO: check if Constraint instance (need ABC module)

        self.constraints.append(constraint)

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """
        """

        observable = np.multiply.reduce(
                [constraint.get(source_coord, telescope) \
                 for constraint in self.constraints])

        return observable

#==============================================================================

class Constraint(object, metaclass=ABCMeta):
    """
    """

    #--------------------------------------------------------------------------
    @abstractmethod
    def __init__(self):
        """
        """

    #--------------------------------------------------------------------------
    @abstractmethod
    def get(self, source_coord, telescope):
        """
        """

#==============================================================================

class ElevationLimit(Constraint):
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit):
        """
        """

        self.limit = limit * u.deg

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """
        """

        altaz = source_coord.transform_to(telescope.frame)
        observable = altaz.alt >= self.limit
        observable = np.array(observable, dtype=float)

        return observable

#==============================================================================

class AirmassLimit(Constraint):
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit, conversion="secz"):
        """
        """

        self.limit = limit

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """
        """

        altaz = source_coord.transform_to(telescope.frame)
        observable = ut.alt_to_airmass(altaz.alt) <= self.limit
        observable = np.array(observable, dtype=float)

        return observable

#==============================================================================

class SunDistance(Constraint):
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit):
        """
        """

        self.limit = limit * u.deg

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """
        """

        source_altaz = source_coord.transform_to(telescope.frame)
        sun_altaz = get_sun(telescope.time).transform_to(telescope.frame)
        separation = source_altaz.separation(sun_altaz)
        observable = separation > self.limit
        observable = np.array(observable, dtype=float)

        return observable

#==============================================================================

class MoonDistance(Constraint):
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit):
        """
        """

        self.limit = limit * u.deg

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """
        """

        source_altaz = source_coord.transform_to(telescope.frame)
        moon_altaz = get_moon(telescope.time).transform_to(telescope.frame)
        separation = source_altaz.separation(moon_altaz)
        observable = separation > self.limit
        observable = np.array(observable, dtype=float)

        return observable

#==============================================================================

class MoonPolarization(Constraint):
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self, limit):
        """
        """

        self.limit = limit * u.deg
        self.limit_lo = (90. - limit) * u.deg
        self.limit_hi = (90. + limit) * u.deg

    #--------------------------------------------------------------------------
    def get(self, source_coord, telescope):
        """
        """

        source_altaz = source_coord.transform_to(telescope.frame)
        moon_altaz = get_moon(telescope.time).transform_to(telescope.frame)
        separation = source_altaz.separation(moon_altaz)
        observable = np.logical_or(
                separation <= self.limit_lo, separation >= self.limit_hi)
        observable = np.array(observable, dtype=float)

        return observable
