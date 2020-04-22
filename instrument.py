# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Classes for instruments.
"""

from abc import ABCMeta, abstractmethod
from astropy.time import TimeDelta
import astropy.units as u

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

class Instrument(object, metaclass=ABCMeta):
    """
    """

    #--------------------------------------------------------------------------
    @abstractmethod
    def __init__(self, readout_time, name=''):
        """
        """

    #--------------------------------------------------------------------------
    @abstractmethod
    def get_obs_time(self, exposure_time, exposure_rep):
        """
        """

#==============================================================================

class InstrumentSimple(Instrument):
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self, readout_time, overhead_time, name=''):
        """
        """

        self.name = name
        self.readout_time = readout_time * u.s
        self.overhead_time = overhead_time * u.s

    #--------------------------------------------------------------------------
    def get_obs_time(self, exposure_time, exposure_rep):
        """
        """

        obs_time = (exposure_time + self.readout_time + self.overhead_time) \
                * exposure_rep
        return obs_time