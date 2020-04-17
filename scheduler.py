# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Scheduler classes.
"""

from abc import ABCMeta, abstractmethod
from astropy.coordinates import get_sun
from astropy.time import Time
import astropy.units as u
import numpy as np

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

class Scheduler(object, metaclass=ABCMeta):
    """
    """

    #--------------------------------------------------------------------------
    @abstractmethod
    def __init__(self):
        """
        """

        self.time_start = None
        self.time_stop = None

    #--------------------------------------------------------------------------
    def _get_time_range(
            self, telescope, twilight, start_time, time_frame):
        """
        """

        # day observations:
        if not twilight:
            # TODO
            raise NotImplementedError()

        # night observations:
        else:
            self.time_start, self.time_stop = telescope.next_sun_set_rise(
                twilight)

    #--------------------------------------------------------------------------
    @abstractmethod
    def run(self, telescope, sources, constraints, duration=None,
            twilight='astronomical', start_time=None, time_frame='utc'):
        """
        """

        # TODO: allow following options:
        # duration: if None, schedule as many days as needed to observe all
        #    sources, otherwise give integer for number of days
        # twilight: if False, schedule full day, otherwise schedule only
        #    observations at night with the following options defining start
        #    and end of the night: 'astronomical' (-18 deg Sun altitude),
        #    'nautical' (-12 deg), 'civic' (-6 deg), 'sunset' (0 deg), or float
        #    to set a specific Sun altitude limit
        # start_time: e.g. 22:00:00;
        #    day observations: start at that time
        #    night observations: start_time is ignored
        # time_frame:
        #    'utc' default
        #    'lst' only works for day observations, sidereal days are used,
        #          start_time is interpreted as sidereal hour

        # program outline:
        # 1. get start time based on sunset, start_time, time_frame
        # 2. get end time if duration
        # 3. Check which targets are visible given the constraints
        # 4. Determine first target
        # 5. Calculate time needed to slew to and observe target
        # 6. Set new time
        # 7. Iterate through steps 3-6 until end time is reached or all sources
        # are scheduled

        # I will need a Schedule class that stores the schedule and can provide
        # info about the schedule, e.g. total time, total slew time, total obs
        # time, total wait time, number of sources, etc...


#==============================================================================

class SimpleScheduler(Scheduler):
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self):
        """
        """

        self.time_start = None
        self.time_stop = None

    #--------------------------------------------------------------------------
    def run(self, telescope, sources, constraints, duration=None,
            twilight='astronomical', start_time=None, time_frame='utc'):
        """
        """

        self._get_time_range(
            telescope, twilight, start_time, time_frame)

#==============================================================================

