# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Scheduler classes.
"""

from abc import ABCMeta, abstractmethod

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

    #--------------------------------------------------------------------------
    @abstractmethod
    def run(self, telescope, sources, constraints, sunset='astronomical',
            duration='day', start_time=None, time_frame='utc'):
        """
        """

        # TODO: allow following options:
        # sunset: 'astronomical', 'nautical', 'civic', float, False
        # duration: 'day' (schedule one day/night), None/False/0 (as many
        # days/nights as necessary to schedule all sources)
        # start_time: e.g. 22:00:00; is overwritten when sunset is defined.
        # time_frame: 'utc', 'local','lst'; defines which frame is used for
        # start_time

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

#==============================================================================

