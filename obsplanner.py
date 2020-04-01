# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Base module of the obsplanner. Provides Obsplanner class.
"""

from astropy.coordinates import SkyCoord
import numpy as np

from constraints import Constraints
from constraints import AirmassLimit, SunDistance, MoonDistance, MoonPolarization
from scheduler import SimpleScheduler
from sources import Sources
from telescope import TelescopeEq, SlewModelLinear

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

class Obsplanner(object):
    """
    """

    #--------------------------------------------------------------------------
    def __init__(self):
        """Create Obsplanner instance.

        Returns
        -----
        None
        """

        self.telescope = None
        self.sources = None
        self.constraints = Constraints()
        self.scheduler = None

    #--------------------------------------------------------------------------
    def set_telescope(self, telescope):
        """Set the telescope for observations.

        Parameters
        -----
        telescope : Telescope instance
            Defines the motion of the telescope.
        """

        # TODO: check if Telescope instance (need ABC module)

        ready = telescope.is_set()
        if ready:
            self.telescope = telescope
            print('Obsplanner: Telescope added.')
        else:
            print('Obsplanner: WARNING: Cannot add Telescope instance.' \
                  'Telescope is not completely set up with slew models.')

    #--------------------------------------------------------------------------
    def set_sources(self, sources):
        """Set the source list for observations.

        Parameters
        -----
        sources : Sources instance
            Sources that should be scheduled.
        """

        # TODO: check if Sources instance (need ABC module)

        self.sources = sources
        print('Obsplanner: Sources added.')

    #--------------------------------------------------------------------------
    def set_scheduler(self, scheduler):
        """Set the scheduler that is used for scheduling the sources.

        Parameters
        -----
        scheduler : Scheduler instance
            Defines the routines how sources are scheduled.
        """

        # TODO: check if Scheduler instance (need ABC module)

        self.scheduler = scheduler

    #--------------------------------------------------------------------------
    def init_time(self, time):
        """Set initial time.

        Parameters
        -----
        time : astropy.Time
            Date and time in UTC. Accepts any input format that astropy. Time
            is accepting.

        Returns
        -----
        None
        """

        if self.telescope is None:
            print('Obsplanner: WARNING: No time set. Set telescope first.')
        else:
            self.telescope.set_time(time)
            print('Obsplanner: Initial date+time:',
                  self.telescope.time)

    #--------------------------------------------------------------------------
    def init_pos(self, coord='zenith'):
        """Set initial telescope position.

        Parameters
        -----
        coord : astropy.Coord or str, default='zenith'
            Specify inital coordinates or set telescope to zenith.

        Returns
        -----
        None
        """

        if self.telescope is None:
            print('Obsplanner: WARNING: No telescope position set. ' \
                  'Set telescope first.')
        elif self.telescope.time is None:
            print('Obsplanner: WARNING: No telescope position set.' \
                  'Set initial date/time first.')
        elif coord=='zenith':
            self.telescope.set_to_zenith()
            print('Obsplanner: Telescope inital position set to zenith.')
        else:
            self.telescope.set_pos(coord)
            print('Obsplanner: Telescope inital position set to:')
            print(self.telescope.pos)


    #--------------------------------------------------------------------------
    def add_constraint(self, constraint):
        """Add an observational constraint.

        Parameters
        -----
        constraint : Constraint instance
            Add a constraint that defines whether or not sources are observable
            at a given time.
        """

        # TODO: chack that constraint is a Constraint instance

        self.constraints.add(constraint)

        # TODO: write which type of constraint was added

    #--------------------------------------------------------------------------
    def _is_ready(self):
        """Check if obsplanner is ready to run the scheduler.
        """

        ready = True

        if self.telescope is None:
            print('Obsplanner: WARNING: Cannot start scheduler. No telescope' \
                  ' set.')
            ready = False
        else:
            if self.telescope.time is None:
                print('Obsplanner: WARNING: Cannot start scheduler. No ' \
                      'initial date/time set.')
                ready = False
            elif self.telescope.pos is None:
                self.init_pos()

        if self.sources is None:
            print('Obsplanner: WARNING: Cannot start scheduler. No sources' \
                  ' set.')
            ready = False

        if self.scheduler is None:
            print('Obsplanner: WARNING: Cannot start scheduler. No scheduler' \
                  ' set.')
            ready = False

        return ready

    #--------------------------------------------------------------------------
    def run_scheduler(self):
        """Run the scheduler.
        """

        if not self._is_ready():
            return False

        print('\nStarting scheduler..')
        self.scheduler.run(
                self.telescope, self.sources, self.constraints,
                sunset='astronomical', duration='day', start_time=None,
                time_frame='utc')



#==============================================================================
# MAIN
#==============================================================================

if __name__ == "__main__":

    obsplanner = Obsplanner()

    telescope = TelescopeEq(
            '24:53:57 deg', '35:12:43 deg', 1750, name='Skinakas')
    slewmodel_ra = SlewModelLinear(2., 0.5, 2.1, 0.6)
    slewmodel_dec = SlewModelLinear(1.9, 0.4, 1.8, 0.5)
    slewmodel_dome = SlewModelLinear(20., 0.1)
    telescope.set_slew_model('ra', slewmodel_ra)
    telescope.set_slew_model('dec', slewmodel_dec)
    #telescope.set_slew_model('dome', slewmodel_dome)
    obsplanner.set_telescope(telescope)
    obsplanner.init_time('2019-11-01 00:00:00')
    #obsplanner.init_pos()
    coord = SkyCoord('75d12m14.1035s', '24d53m57s')
    #obsplanner.init_pos(coord)

    dtype = [('name', 'S30'), ('ra', float), ('dec', float), ('expt', float),
             ('expn', int)]
    targets = np.loadtxt(
            'sourcelists/targets.csv', dtype=dtype, skiprows=1, delimiter=',',
            usecols=range(5))
    sources = Sources(targets['name'], targets['ra'], targets['dec'])
    obsplanner.set_sources(sources)

    airmass_limit = AirmassLimit(2.)
    sun_distance = SunDistance(20.)
    moon_distance = MoonDistance(20.)
    moon_polarization = MoonPolarization(10.)

    obsplanner.add_constraint(airmass_limit)
    obsplanner.add_constraint(sun_distance)
    obsplanner.add_constraint(moon_distance)
    obsplanner.add_constraint(moon_polarization)

    scheduler = SimpleScheduler()
    obsplanner.set_scheduler(scheduler)
    obsplanner.run_scheduler()



    #telescope.set_time('2019-11-01 00:00:00')
    #telescope.set_to_zenith()
    # TODOs:
    # 1. I do not want to set the time for the telescope but for the obsplanner
    # 2. If not start position is applied, the telescope should be set to
    #    zenith automatically, before running the scheduler
    # 3. I want it flexible to set first the start position then the time or
    #    vice versa.
    # When the scheduler sets the start time to night start, the initial
    # position needs to be updated.

    # should I add the concept of preferences? similar to Constraints, but
    # returning values between 0 and 1