# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Observation block classes.
"""

from astropy import units as u
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
    """
    """

    #--------------------------------------------------------------------------
    def __init__(
            self, source_id, source_name, source_coord, exp_time, exp_rep,
            time_start, time_slew, time_obs, time_tot):
        """


        Parameters
        ----------
        source_id : TYPE
            DESCRIPTION.
        source_name : TYPE
            DESCRIPTION.
        source_coord : TYPE
            DESCRIPTION.
        exp_time : TYPE
            DESCRIPTION.
        exp_rep : TYPE
            DESCRIPTION.
        time_start : TYPE
            DESCRIPTION.
        time_slew : TYPE
            DESCRIPTION.
        time_obs : TYPE
            DESCRIPTION.
        time_tot : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

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
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self, block_id, time_start, time_stop):
        """

        Parameters
        ----------
        time_start : TYPE
            DESCRIPTION.
        time_stop : TYPE
            DESCRIPTION.

        Returns
        -------
        None.
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
    def add_observation(
            self, source_id, source_name, source_coord, source_exp, source_rep,
            time_start, time_slew, time_obs, time_tot):
        """


        Parameters
        ----------
        source_id : TYPE
            DESCRIPTION.
        source_name : TYPE
            DESCRIPTION.
        source_coord : TYPE
            DESCRIPTION.
        source_exp : TYPE
            DESCRIPTION.
        source_rep : TYPE
            DESCRIPTION.
        time_now : TYPE
            DESCRIPTION.
        time_slew : TYPE
            DESCRIPTION.
        time_obs : TYPE
            DESCRIPTION.
        time_tot : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

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
        """


        Returns
        -------
        info : TYPE
            DESCRIPTION.

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

#==============================================================================

class Schedule(object):
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self):
        """

        Returns
        -------
        None.

        """

        self.blocks = []
        self.last_block = False
        self.n_blocks = 0
        self.n_obs = 0

    #--------------------------------------------------------------------------
    def add_obsblock(self, time_start, time_stop):
        """


        Parameters
        ----------
        block : TYPE
            DESCRIPTION.
        time_start : TYPE
            DESCRIPTION.
        time_stop : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        self.blocks.append(ObsBlock(self.n_blocks, time_start, time_stop))
        self.last_block = self.blocks[-1]
        self.n_blocks += 1

    #--------------------------------------------------------------------------
    def add_observation(
            self, source_id, source_name, source_coord, source_exp, source_rep,
            time_start, time_slew, time_obs, time_tot):
        """


        Parameters
        ----------
        source_id : TYPE
            DESCRIPTION.
        source_name : TYPE
            DESCRIPTION.
        source_coord : TYPE
            DESCRIPTION.
        source_exp : TYPE
            DESCRIPTION.
        source_rep : TYPE
            DESCRIPTION.
        time_start : TYPE
            DESCRIPTION.
        time_slew : TYPE
            DESCRIPTION.
        time_obs : TYPE
            DESCRIPTION.
        time_tot : TYPE
            DESCRIPTION.

        Raises
        ------
        AttributeError
            DESCRIPTION.

        Returns
        -------
        None.

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
        """


        Yields
        ------
        TYPE
            DESCRIPTION.

        """

        if self.n_blocks == 0:
            return False

        for block in self.blocks:
            yield block
