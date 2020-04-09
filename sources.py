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

        self.name = np.array(name)

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
        self.exposure_time = TimeDelta(exposure_time, format='sec')
        self.exposure_rep = exposure_rep
        self.active = np.ones(self.coord.size, dtype=bool)
        self.scheduled = np.ones(self.coord.size, dtype=bool)
        self.size = self.name.size

    #--------------------------------------------------------------------------
    def __len__(self):
        return self.size
        
    

