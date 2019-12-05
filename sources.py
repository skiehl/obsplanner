# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Providing the Sources class.
"""

import numpy as np
from astropy.coordinates import Angle, SkyCoord
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
            self, name, ra, dec, ra_unit='deg', dec_unit='deg'):
        """Create Sources instance.

        Parameters
        -----
        name : list of str
            Source names.
        ra : np.ndarray
            Array of right ascension values.
        dec : np.ndarray
            Array of declination values.
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

