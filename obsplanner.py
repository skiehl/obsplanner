# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Base module of the obsplanner. Provides Obsplanner class.
"""

from astropy.coordinates import SkyCoord
import numpy as np

from constraints import Constraint, Constraints
from constraints import AirmassLimit, ElevationLimit, SunDistance, MoonDistance, MoonPolarization
from instrument import Instrument, InstrumentSimple
from scheduler import Scheduler, SimpleScheduler
from sources import Sources
from telescope import Telescope, TelescopeEq, SlewModelLinear

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
    """Observation planner.

    Combines all components needed for observation planning in order to run the
    planner.
    """

    #--------------------------------------------------------------------------
    def __init__(self):
        """Create Obsplanner instance.

        Returns
        -----
        None
        """

        self.instrument = None
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

        # check if Telescope instance:
        if not isinstance(telescope, Telescope):
            raise TypeError('Unsupported type: {0}'.format(type(telescope)))

        # check if telescope is readily set up:
        ready = telescope.is_set()
        if ready:
            self.telescope = telescope
            if telescope.name is None:
                print('Obsplanner: Telescope added.')
            else:
                print('Obsplanner: Telescope {0:s} added.'.format(
                        telescope.name))
        else:
            print('Obsplanner: WARNING: Cannot add Telescope instance.' \
                  'Telescope is not completely set up with slew models.')

    #--------------------------------------------------------------------------
    def set_instrument(self, instrument):
        """Set the instrument for observations.

        Parameters
        -----
        instrument : Insturment instance
            Defines the time needed for each observation.
        """

        # check if Instrument instance:
        if not isinstance(instrument, Instrument):
            raise TypeError('Unsupported type: {0}'.format(type(instrument)))

        self.instrument = instrument
        if instrument.name is None:
            print('Obsplanner: Instrument added.')
        else:
            print('Obsplanner: Instrument {0:s} added.'.format(
                    instrument.name))

    #--------------------------------------------------------------------------
    def set_sources(self, sources):
        """Set the source list for observations.

        Parameters
        -----
        sources : Sources instance
            Sources that should be scheduled.
        """

        # check if Sources instance:
        if not isinstance(sources, Sources):
            raise TypeError('Unsupported type: {0}'.format(type(sources)))

        self.sources = sources
        print('Obsplanner: {0:d} sources added.'.format(sources.size))

    #--------------------------------------------------------------------------
    def set_scheduler(self, scheduler):
        """Set the scheduler that is used for scheduling the sources.

        Parameters
        -----
        scheduler : Scheduler instance
            Defines the routines how sources are scheduled.
        """

        # check if Scheduler instance:
        if not isinstance(scheduler, Scheduler):
            raise TypeError('Unsupported type: {0}'.format(type(scheduler)))

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
            print('Obsplanner: Initial date+time:', time)
            self.telescope.set_time(time)

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

        self.constraints.add(constraint)
        print('Obsplanner: constraint added: {0:s}'.format(
                constraint.__str__()))

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

        if self.instrument is None:
            print('Obsplanner: WARNING: Cannot start scheduler. No ' \
                  'instrument set.')
            ready = False

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
    def run_scheduler(self, duration=None):
        """Run the scheduler.
        """

        if not self._is_ready():
            return False

        print('\nStarting scheduler..')
        self.scheduler.run(
                self.telescope, self.instrument, self.sources,
                self.constraints, duration=duration, twilight='astronomical',
                start_time=None, time_frame='utc')

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
    obsplanner.init_time('2019-06-01 00:00:00')
    #obsplanner.init_pos()
    coord = SkyCoord('75d12m14.1035s', '24d53m57s')
    #obsplanner.init_pos(coord)

    instrument = InstrumentSimple(100., 0., name='WALOP')
    obsplanner.set_instrument(instrument)

    dtype = [('name', 'S30'), ('ra', float), ('dec', float), ('expt', float),
             ('expn', int)]
    targets = np.loadtxt(
            'unittests/sourcelists/targets_fewer.csv', dtype=dtype, skiprows=1,
            delimiter=',', usecols=range(5))
    sources = Sources(
            targets['name'], targets['ra'], targets['dec'], targets['expt'],
            targets['expn'])
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
    obsplanner.run_scheduler(duration=None)



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