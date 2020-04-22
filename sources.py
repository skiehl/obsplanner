# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Providing the Sources class.
"""

import numpy as np
from astropy.coordinates import Angle, SkyCoord
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

class Sources:
    """List of sources.
    """

    #--------------------------------------------------------------------------
    def __init__(
            self, name, ra, dec, exposure_time, exposure_rep, ra_unit='deg',
            dec_unit='deg'):
        """Create Sources instance.

        Parameters
        -----
        name : list of str
            Source names.
        ra : np.ndarray
            Array of right ascension values.
        dec : np.ndarray
            Array of declination values.
        exposure_time : np.ndarray
            Exposure times for each source in seconds.
        exposure_rep : np.ndarray
            Number of exposure repetitions for each source.
        ra_unit : str, default='deg'
            Unit right ascensions are provided in.
        ra_unit : str, default='deg'
            Unit declinations are provided in.

        Returns
        -----
        None

        Notes
        -----
        Currently only degrees i.e. 'deg' is provided as option. Further
        options may be implemented later.
        """

        self.id = np.arange(len(name))
        #self.name = np.array(name)
        self.name = name

        # store right ascension:
        if ra_unit == 'deg':
            ra = Angle(ra * u.deg)
        else:
            raise ValueError("ra_unit needs to be 'deg'.")

        # store declination:
        if dec_unit == 'deg':
            dec = Angle(dec * u.deg)
        else:
            raise ValueError("dec_unit needs to be 'deg'.")

        self.coord = SkyCoord(ra, dec)
        self.exposure_time = exposure_time * u.s
        self.exposure_rep = exposure_rep
        self.active = np.ones(self.coord.size, dtype=bool)
        self.scheduled = np.zeros(self.coord.size, dtype=bool)
        self.observable_now = np.ones(self.coord.size, dtype=bool)
        self.observable_general = np.zeros(self.coord.size, dtype=bool)
        self.size = self.name.size

        print('Sources: Source list created with {0:d} sources.'.format(
                self.size))

    #--------------------------------------------------------------------------
    def __len__(self):
        return self.size

    #--------------------------------------------------------------------------
    def set_observable(self, source_id, observable, verbose=True):
        """Set which sources are currently observable.

        Parameters
        ----------
        source_id : np.ndarray (int)
            IDs of the sources that are currently observable.
        observable : numpy.ndarray (bool)
            Each entry corresponds to one source.
            True, if observable. False, otherwise.
        verbose : bool, default=True
            If True, print out information.

        Returns
        -------
        None.
        """

        self.observable_now = np.zeros(self.size, dtype=bool)
        self.observable_now[source_id] = observable
        self.observable_general = np.logical_or(
                self.observable_general, self.observable_now)

        if verbose:
            print('Sources: {0:d} sources set not observable.'.format(
                    self.size-np.sum(observable)))

    #--------------------------------------------------------------------------
    def set_scheduled(self, source_id):
        """Set source(s) as scheduled.

        Parameters
        ----------
        source_id : np.ndarray (int) or int
            ID(s) of the source(s) that has/have been scheduled.

        Returns
        -------
        None.
        """

        self.scheduled[source_id] = True

    #--------------------------------------------------------------------------
    def get_source(self, source_id):
        """Return source parameters.

        Parameters
        ----------
        source_id : int
            ID of the source.

        Returns
        -------
        str
            Source name.
        SkyCoord
            Source coordinates.
        float
            Exposure time
        int
            Exposure repetitions
        """

        return (
            self.name[source_id], self.coord[source_id],
            self.exposure_time[source_id], self.exposure_rep[source_id])

    #--------------------------------------------------------------------------
    def _select(
            self, active, scheduled, observable_now, observable_general,
            verbose):
        """Select sources according to certain selection criteria.

        Parameters
        ----------
        active : bool or None
            If True, select only active sources.
            If False, select only inactive sources.
            If None, select both.
        scheduled : bool or None
            If True, select only sources that have been scheduled already.
            If False, select only sources that are not yet scheduled.
            If None, select both.
        observable_now : bool or None
            If True, select only sources currently observable.
            If False, select only sources currently not observable.
            If None, select both.
        observable_general : bool or None
            If True, select only sources that were observable at some time.
            If False, select only sources that were not observable at any time.
            If None, select both.
        verbose : bool, default=True
            If True, print out information.

        Returns
        -------
        np.ndarray
            Array of the same size as the number sources with True entries for
            those selected.

        Notes
        -------
        Used by self.get_sources() and self.count_sources().
        """



        info = []
        if active is True:
            sel = self.active
            info.append('active')
        elif active is False:
            sel = ~self.active
            info.append('inactive')
        else:
            sel = np.ones(self.size, dtype=bool)

        if scheduled is True:
            sel = np.logical_and(sel, self.scheduled)
            info.append('already scheduled')
        elif scheduled is False:
            sel = np.logical_and(sel, ~self.scheduled)
            info.append('not yet scheduled')
        else:
            pass

        if observable_now is True:
            sel = np.logical_and(sel, self.observable_now)
            info.append('currently observable')
        elif observable_now is False:
            sel = np.logical_and(sel, ~self.observable_now)
            info.append('currently not observable')
        else:
            pass

        if observable_general is True:
            sel = np.logical_and(sel, self.observable_general)
            info.append('generally observable')
        elif observable_general is False:
            sel = np.logical_and(sel, ~self.observable_general)
            info.append('generally not observable')
        else:
            pass

        if verbose:
            info_text = 'Sources: {0:d}'.format(np.sum(sel))
            for el in info:
                info_text = '{0:s} {1:s},'.format(info_text, el)
            if len(info):
                info_text = info_text[:-1]
            info_text = '{0:s} sources selected.'.format(info_text)
            print(info_text)

        return sel

    #--------------------------------------------------------------------------
    def get_sources(
            self, active=None, scheduled=None, observable_now=None,
            observable_general=None, verbose=True):
        """Return source IDs, names, coordinates, exposure times and
        repetiotions.

        Parameters
        ----------
        active : bool or None, default=None
            If True, return only active sources.
            If False, return only inactive sources.
            If None, return both.
        scheduled : bool or None, default=None
            If True, return only sources that have been scheduled already.
            If False, return only sources that are not yet scheduled.
            If None, return both.
        observable_now : bool or None, default=None
            If True, return only sources currently observable.
            If False, return only sources currently not observable.
            If None, return both.
        observable_general : bool or None, default=None
            If True, return only sources that were observable at some time.
            If False, return only sources that were not observable at any time.
            If None, return both.
        verbose : bool, default=True
            If True, print out information.

        Returns
        -------
        id : np.ndarray
            Source IDs of the selected sources.
        coord :astropy.coordinates.sky_coordinate.SkyCoord
            Coordinates of selected sources.
        name : numpy.ndarray
            Names of selected sources.
        exposure_time : numpy.ndarray
            Exposure times for selected sources.
        exposure_rep : numpy.ndarray
            Exposure repetitions for selected sources.

        Notes
        -------
        active: only active sources are supposed to get scheduled.
        scheduled: scheduled means the source has already been scheduled.
        observable_now: sources that are currently observable given certain
        observational constraints.
        observable_general: sources that were observable at least at one moment
        during a scheduled time period given certain observational constraints.
        """

        sel = self._select(
            active, scheduled, observable_now, observable_general, verbose)

        if np.sum(sel) == 0:
            return False

        return (
            self.id[sel], self.coord[sel], self.name[sel],
            self.exposure_time[sel], self.exposure_rep[sel])

    #--------------------------------------------------------------------------
    def count_sources(
            self, active=None, scheduled=None, observable_now=None,
            observable_general=None):
        """Count sources with given criteria.

        Parameters
        ----------
        active : bool or None, default=None
            If True, count only active sources.
            If False, count only inactive sources.
            If None, count both.
        scheduled : bool or None, default=None
            If True, count only sources that have been scheduled already.
            If False, count only sources that are not yet scheduled.
            If None, count both.
        observable_now : bool or None, default=None
            If True, count only sources currently observable.
            If False, count only sources currently not observable.
            If None, count both.
        observable_general : bool or None, default=None
            If True, count only sources that were observable at some time.
            If False, count only sources that were not observable at any time.
            If None, count both.

        Returns
        -------
        int
            Number of selected sources.

        Notes
        -------
        active: only active sources are supposed to get scheduled.
        scheduled: scheduled means the source has already been scheduled.
        observable_now: sources that are currently observable given certain
        observational constraints.
        observable_general: sources that were observable at least at one moment
        during a scheduled time period given certain observational constraints.
        """

        sel = self._select(
            active, scheduled, observable_now, observable_general, False)
        n = np.sum(sel)

        return n


