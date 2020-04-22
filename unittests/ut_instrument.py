# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Module description.
"""

import sys
sys.path.append('../')
from astropy import units as u

import instrument as i

__author__ = "Sebastian Kiehlmann"
__credits__ = ["Sebastian Kiehlmann"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sebastian Kiehlmann"
__email__ = "skiehlmann@mail.de"
__status__ = "Production"

#==============================================================================
# MAIN
#==============================================================================

print('\nTEST instrument.InstrumentSimple:')
instrument = i.InstrumentSimple(100., 0., name='RoboPol')

#------------------------------------------------------------------------------
print('\nTEST instrument.InstrumentSimple.get_obs_time():')
exposure_time = 20. * u.s
exposure_rep = 5
print(instrument.get_obs_time(exposure_time, exposure_rep))
