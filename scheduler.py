# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Scheduler classes.
"""

from abc import ABCMeta, abstractmethod
from astropy.coordinates import get_sun
from astropy.time import Time
import astropy.units as u
import numpy as np

import obsblock as o

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
        self.lst_start = None
        self.lst_stop = None

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
            self.lst_start = self.time_start.sidereal_time('apparent')
            self.lst_stop = self.time_stop.sidereal_time('apparent')

    #--------------------------------------------------------------------------
    def _observable_sources(
            self, telescope, sources, constraints, interval=5.):
        """Determine the observable sources at a given time.

        Parameters
        -------
        telescope : Telescope instance
            Defines the telescope location and current time.
        sources : Sources instance
            Sources that should be scheduled.
        constraints : Constraints instance
            Defines observational constraints.
        interval : float, default=5.
            Time interval in minutes. If no sources are observable at a given
            time, the time is continuously increased by this interval until new
            sources become observable.


        Returns
        -------
        time_now : Time
            The current time of the telescope at which the next sources are
            observable.
        """

        # find observable sources:
        source_id, source_coord, __, __, __ = sources.get_sources(
                active=True, scheduled=False)
        observable = constraints.get(source_coord, telescope)
        sources.set_observable(source_id, observable)
        n = np.sum(observable)

        time_now = telescope.get_time()

        if n == 0:
            print('Scheduler: No sources observable at', time_now)
            print('  Advancing in time..')

        # update time, if no sources are observable:
        while n == 0 and time_now < self.time_stop:
            # advance time:
            time_now = time_now + interval * u.min
            telescope.set_time(time_now, verbose=False)

            # find observable sources:
            source_id, source_coord, __, __, __ = sources.get_sources(
                    active=True, scheduled=False, verbose=False)
            observable = constraints.get(source_coord, telescope)
            sources.set_observable(source_id, observable, verbose=False)
            n = np.sum(observable)

        if n > 0:
            print('  Continue at', time_now)
            return time_now
        else:
            return False

    #--------------------------------------------------------------------------
    @abstractmethod
    def run(self, telescope, instrument, sources, constraints, duration=None,
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
        self.schedule = None

    #--------------------------------------------------------------------------
    def _schedule_block(
            self, telescope, instrument, sources, constraints, twilight,
            start_time, time_frame):
        """
        """

        # set telescope to start time:
        self._get_time_range(
            telescope, twilight, start_time, time_frame)
        telescope.set_time(self.time_start)

        # create observation block:
        self.schedule.add_obsblock(self.time_start, self.time_stop)

        # count observable sources:
        n_iter = sources.count_sources(active=True, scheduled=False)

        # iterate for number of sources:
        for __ in range(n_iter):
            time_now = self._observable_sources(
                telescope, sources, constraints, interval=5.)
            if time_now is False:
                print('Scheduler: end of observation window reached.')
                break

            # get slew time:
            sources_id, sources_coord, sources_name, sources_exp, sources_rep \
                    = sources.get_sources(
                        active=True, scheduled=False, observable_now=True)
            slew_time = telescope.get_slew_time(sources_coord)

            # select shortest slew time:
            i = np.argmin(slew_time)
            time_slew = slew_time[i]
            time_obs = instrument.get_obs_time(sources_exp[i], sources_rep[i])
            time_tot = time_slew + time_obs
            time_new = time_now + time_tot

            # stop if new time exceeds stop time:
            if time_new > self.time_stop:
                print('Scheduler: end of observation window reached.')
                break

            # add source to schedule:
            if isinstance(sources_name[i], bytes):
                source_name = sources_name[i].decode('UTF-8')
            else:
                source_name = sources_name[i]
            self.schedule.add_observation(
                i, source_name, sources_coord[i],
                sources_exp[i], sources_rep[i], time_now, time_slew, time_obs,
                time_tot)

            # move telescope, set new time:
            telescope.set_pos(sources_coord[i])
            telescope.set_time(time_new)
            sources.set_scheduled(sources_id[i])

        block = self.schedule.last_block
        if block.n_obs == 0:
            return False

        print(block.summary())
        return True

    #--------------------------------------------------------------------------
    def run(self, telescope, instrument, sources, constraints, duration=None,
            twilight='astronomical', start_time=None, time_frame='utc'):
        """
        """

        # TODO: can I move code dublications to other methods?

        # check if schedule already exists:
        if self.schedule is not None:
            user_in = input("Overwrite existing schedule? y/n")
            if user_in in ['y', 'yes', 'Y', 'Yes', 'Make it so!']:
                pass
            else:
                return False

        self.schedule = o.Schedule()

        # schedule specific number of days:
        if duration:
            print('Scheduler: schedule {0:d} days'.format(duration))

            # iterate through days:
            for i in range(duration):
                print('Scheduler: schedule day {0:d}'.format(i+1))
                done = self._schedule_block(
                    telescope, instrument, sources, constraints, twilight,
                    start_time, time_frame)

                if not done:
                    print('Scheduler: no sources scheduled for this block. ' \
                          'Abort.')
                    break

        # schedule all sources:
        else:
            i = 0
            while sources.count_sources(scheduled=False) > 0:
                print('Scheduler: schedule day {0:d}'.format(i+1))
                i += 1
                done = self._schedule_block(
                    telescope, instrument, sources, constraints, twilight,
                    start_time, time_frame)

                if not done:
                    print('Scheduler: no sources scheduled for this block. ' \
                          'Abort.')
                    break


#==============================================================================

