# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Scheduler classes.
"""

from abc import ABCMeta, abstractmethod
import astropy.units as u
import numpy as np
import textwrap

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

class Observation(object):
    """Observation of one source.
    """

    #--------------------------------------------------------------------------
    def __init__(
            self, source_id, source_name, source_coord, exp_time, exp_rep,
            time_start, time_slew, time_obs, time_tot):
        """Observation of one source.

        Parameters
        ----------
        source_id : int
            Source ID.
        source_name : str
            Source name.
        source_coord : astropy.coordinates.SkyCoord
            Source coordinates.
        exp_time : float
            Duration of a single exposure.
        exp_rep : int
            Number of exposure repetitions.
        time_start : astropy.time.Time
            Time when moving to the source for an observation starts.
        time_slew : astropy.units.quantity.Quantity
            Slew time in seconds.
        time_obs : astropy.units.quantity.Quantity
            Observation time in seconds.
        time_tot : astropy.units.quantity.Quantity
            Slew and observation time in seconds.

        Returns
        -------
        None
        """

        self.source_id = source_id
        self.source_name = source_name
        self.source_coord = source_coord
        self.exp_time = exp_time
        self.exp_rep = exp_rep
        self.time_start = time_start
        self.time_slew = time_slew
        self.time_obs = time_obs
        self.time_tot = time_tot

#==============================================================================

class ObsBlock(object):
    """Observation block.
    """

    #--------------------------------------------------------------------------
    def __init__(self, block_id, time_start, time_stop):
        """Observation block.

        Parameters
        ----------
        block_id : int
            Observation block ID.
        time_start : astropy.time.Time
            Time when the observation block starts.
        time_stop : astropy.time.Time
            Time when the observation block ends.

        Returns
        -------
        None
        """

        self.id = block_id
        self.time_start = time_start
        self.time_stop = time_stop
        self.duration = (time_stop - time_start).to_value('sec') * u.s
        self._time_current = time_start
        self.time_slew = 0 * u.s
        self.time_obs = 0 * u.s
        self.time_idle = 0 * u.s
        self.observations = []
        self.n_obs = 0

    #--------------------------------------------------------------------------
    def __len__(self):
        """Number of observations.

        Returns
        -------
        None
        """

        return self.n_obs

    #--------------------------------------------------------------------------
    def __str__(self):
        """Information about the observation block.

        Returns
        -------
        str
            Printed information about the observation block.
        """

        return self.summary()

    #--------------------------------------------------------------------------
    def add_observation(
            self, source_id, source_name, source_coord, source_exp, source_rep,
            time_start, time_slew, time_obs, time_tot):
        """Add new observation to the observation block.

        Parameters
        ----------
        source_id : int
            Source ID.
        source_name : str
            Source name.
        source_coord : astropy.coordinates.SkyCoord
            Source coordinates.
        exp_time : float
            Duration of a single exposure.
        exp_rep : int
            Number of exposure repetitions.
        time_start : astropy.time.Time
            Time when moving to the source for an observation starts.
        time_slew : astropy.units.quantity.Quantity
            Slew time in seconds.
        time_obs : astropy.units.quantity.Quantity
            Observation time in seconds.
        time_tot : astropy.units.quantity.Quantity
            Slew and observation time in seconds.

        Returns
        -------
        None
        """

        # create observation:
        observation = Observation(
            source_id, source_name, source_coord, source_exp,
            source_rep, time_start, time_slew, time_obs, time_tot)
        self.observations.append(observation)
        self.n_obs += 1

        # update slew, obs, and idle time
        self.time_slew += time_slew
        self.time_obs += time_obs
        self.time_idle = self.duration - self.time_slew - self.time_obs

        print("Schedule: observation added: {0:s}".format(source_name))

    #--------------------------------------------------------------------------
    def summary(self):
        """Information about the observation block.

        Returns
        -------
        info : str
            Printed information about the observation block.
        """

        slew_ratio = self.time_slew / self.duration * 100.
        obs_ratio = self.time_obs / self.duration * 100.
        idle_ratio = self.time_idle / self.duration * 100.

        info = textwrap.dedent("""
            Observing block {0:d} summary:
            Scheduled sources: {1:d}
            Start time:        {2:s} UTC
            Stop time:         {3:s} UTC
            Duration:          {4:10.2f}
            Slew time:         {5:10.2f} {8:6.1f} %
            Observing time:    {6:10.2f} {9:6.1f} %
            Idle time:         {7:10.2f} {10:6.1f} %
            """.format(
                self.id, self.n_obs, self.time_start.iso, self.time_stop.iso,
                self.duration, self.time_slew, self.time_obs, self.time_idle,
                slew_ratio, obs_ratio, idle_ratio)
            )

        return info

    #--------------------------------------------------------------------------
    def get_observations(self):
        """Iterator over stored observations.

        Yields
        ------
        obsplanner.scheduler.Observation
            Observations.
        """

        for observation in self.observations:
            yield observation

#==============================================================================

class Schedule(object):
    """Observation schedule.
    """

    #--------------------------------------------------------------------------
    def __init__(self):
        """Observation schedule.

        Returns
        -------
        None
        """

        self.blocks = []
        self.last_block = False
        self.n_blocks = 0
        self.n_obs = 0

    #--------------------------------------------------------------------------
    def add_obsblock(self, time_start, time_stop):
        """Add observation block to schedule.


        Parameters
        ----------
        time_start : astropy.time.Time
            Time when the observation block starts.
        time_stop : astropy.time.Time
            Time when the observation block ends.

        Returns
        -------
        None
        """

        self.blocks.append(ObsBlock(self.n_blocks, time_start, time_stop))
        self.last_block = self.blocks[-1]
        self.n_blocks += 1

    #--------------------------------------------------------------------------
    def add_observation(
            self, source_id, source_name, source_coord, source_exp, source_rep,
            time_start, time_slew, time_obs, time_tot):
        """Add new observation to the observation block.

        Parameters
        ----------
        source_id : int
            Source ID.
        source_name : str
            Source name.
        source_coord : astropy.coordinates.SkyCoord
            Source coordinates.
        exp_time : float
            Duration of a single exposure.
        exp_rep : int
            Number of exposure repetitions.
        time_start : astropy.time.Time
            Time when moving to the source for an observation starts.
        time_slew : astropy.units.quantity.Quantity
            Slew time in seconds.
        time_obs : astropy.units.quantity.Quantity
            Observation time in seconds.
        time_tot : astropy.units.quantity.Quantity
            Slew and observation time in seconds.

        Returns
        -------
        None

        Raises
        ------
        AttributeError
            When no observing block has been set yet.
        """

        if self.n_blocks == 0:
            raise AttributeError(
                "Cannot add observations before adding at least one observing"\
                " block.")

        self.blocks[-1].add_observation(
                source_id, source_name, source_coord, source_exp, source_rep,
                time_start, time_slew, time_obs, time_tot)
        self.n_obs += 1

    #--------------------------------------------------------------------------
    def get_blocks(self):
        """Iterator over observing blocks.


        Yields
        ------
        obsplanner.scheduler.ObsBlock
            Observing block.
        """

        for block in self.blocks:
            yield block

#==============================================================================

class Scheduler(object, metaclass=ABCMeta):
    """Scheduler. Finds an optimized order of target sources.

    Notes
    ------
    This is an abstract base class. Different types of schedulers are
    implemented based on this parent class.
    """

    #--------------------------------------------------------------------------
    @abstractmethod
    def __init__(self):
        """Scheduler. Finds an optimized order of target sources.

        Returns
        ------
        None
        """

        self.time_start = None
        self.time_stop = None
        self.lst_start = None
        self.lst_stop = None
        self.schedule = None

    #--------------------------------------------------------------------------
    def _check_schedule(self):
        """Check if schedule exists.

        Returns
        -------
        bool
            True, if the scheduler should proceed to run and create a new
            schedule.
            False, if the scheduler should abort and keep the current schedule.
        """

        if self.schedule is not None:
            user_in = input("Overwrite existing schedule? y/n")
            if user_in in ['y', 'yes', 'Y', 'Yes', 'Make it so!']:
                pass
            else:
                print('Scheduler: aborted!')
                return False

        return True

    #--------------------------------------------------------------------------
    def _get_time_range(
            self, telescope, twilight, start_time, time_frame):
        """Calculate time range for the next observing block.

        Parameters
        ----------
        telescope : obsplanner.telescope.Telescope
            Defines the position of the observatory.
        twilight : str or float
            Select the Sun elevation at which observations should start/end.
            Choose from 'astronomical' (-18 deg), 'nautical' (-12 deg),
            'civil' (-6 deg), 'sunset' (0 deg). Or use float to set Sun
            elevation (in degrees).
        start_time : TYPE
            DESCRIPTION.
        time_frame : TYPE
            DESCRIPTION.

        Raises
        ------
        NotImplementedError
            Raised if twilight is not set or when a time_frame other than
            'utc' is given.

        Returns
        -------
        None
        """

        # TODO: use of time_frame and start_time not implemented yet
        if time_frame != 'utc':
            # TODO
            raise NotImplementedError()

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
        telescope : obsplanner.telescope.Telescope
            Defines the telescope location and current time.
        sources : obsplanner.sources.Sources
            Sources that should be scheduled.
        constraints : obsplanner.constraints.Constraints
            Defines observational constraints.
        interval : float, default=5.
            Time interval in minutes. If no sources are observable at a given
            time, the time is continuously increased by this interval until new
            sources become observable.


        Returns
        -------
        astropy.time.Time
            The current time of the telescope at which the next sources are
            observable. Returns False if no sources are observable.
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
        """TODO
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

class SchedulerNearestNeighbor(Scheduler):
    """Scheduler picking the nearest neighbor.
    """

    #--------------------------------------------------------------------------
    def __init__(self):
        """Scheduler picking the nearest neighbor.

        Returns
        -------
        None
        """

        self.time_start = None
        self.time_stop = None
        self.schedule = None

    #--------------------------------------------------------------------------
    def _schedule_block(
            self, telescope, instrument, sources, constraints, twilight,
            start_time, time_frame):
        """


        Parameters
        ----------
        telescope : obsplanner.telescope.Telescope
            Defines the telescope location and motion.
        instrument : obsplanner.instrument.Instrument
            Defines the time needed for observations.
        sources : obsplanner.sources.Sources
            Sources that should be scheduled.
        constraints : obsplanner.constraints.Constraints
            Defines observational constraints.
        twilight : str or float
            Select the Sun elevation at which observations should start/end.
            Choose from 'astronomical' (-18 deg), 'nautical' (-12 deg),
            'civil' (-6 deg), 'sunset' (0 deg). Or use float to set Sun
            elevation (in degrees).
        start_time : TYPE
            DESCRIPTION.
        time_frame : TYPE
            DESCRIPTION.

        Returns
        -------
        bool
            True, if new observing block has been scheduled. False, otherwise.

        Raises
        ------
        NotImplementedError
            Raised if a time_frame other than 'utc' is given.

        Notes
        ------
        Use of time_frame and start_time is not implemented yet.
        """

        # TODO: use of time_frame and start_time not implemented yet
        if time_frame != 'utc':
            # TODO
            raise NotImplementedError()

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
        """TODO
        """

        # check if schedule exists:
        if not self._check_schedule():
            return False

        self.schedule = Schedule()

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

