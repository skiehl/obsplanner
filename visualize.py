# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""Visualization functions.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
import numpy as np

__author__ = "Sebastian Kiehlmann"
__credits__ = ["Sebastian Kiehlmann"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Sebastian Kiehlmann"
__email__ = "skiehlmann@mail.de"
__status__ = "Production"

#==============================================================================
# FUNCTIONS
#==============================================================================

def ax_add_obs(obs, ax):
    """


    Parameters
    ----------
    obs : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.

    Returns
    -------
    None.
    """

    linewidth = 1
    linestyle = '-'
    edgecolor= 'k'

    # plot slew window:
    x = obs.time_start.to_value('mjd')
    width_slew = obs.time_slew.to_value('d')
    rect = Rectangle(
            (x, 0), width_slew, 1, facecolor='0.7', linestyle='None')
    ax.add_patch(rect)

    # plot observation window:
    x = obs.time_start + obs.time_slew
    x = x.to_value('mjd')
    width_obs = obs.time_obs.to_value('d')
    rect = Rectangle(
            (x, 0), width_obs, 1, facecolor='0.6', linestyle='None')
    ax.add_patch(rect)

    # plot complete window:
    x = obs.time_start.to_value('mjd')
    rect = Rectangle(
            (x, 0), width_slew+width_obs, 1, fill=False,
            linewidth=linewidth, linestyle=linestyle, edgecolor=edgecolor)
    ax.add_patch(rect)

#------------------------------------------------------------------------------
def plot_block(block, utc_offset=0):
    """


    Returns
    -------
    None.
    """

    linewidth = 1
    linestyle = '-'
    edgecolor= 'k'
    y = 1.2

    gridspec_kw = {'width_ratios': (5,1), 'wspace': 0.2}
    fig, ax = plt.subplots(1, 2, figsize=(16,9), gridspec_kw=gridspec_kw)

    # plot block:
    x = block.time_start.to_value('mjd')
    width = block.duration.to_value('d')
    rect = Rectangle(
            (x, 0), width, y, fill=False, linewidth=linewidth,
            linestyle=linestyle, edgecolor=edgecolor)
    ax[0].add_patch(rect)

    for i, obs in enumerate(block.get_observations()):
        ax_add_obs(obs, ax[0])
    xmin = block.time_start.to_value('mjd') - 1. / 24
    xmax = block.time_stop.to_value('mjd') + 1. / 24
    ax[0].set_xlim(xmin, xmax)
    ax[0].xaxis.set_major_locator(MultipleLocator(1./12))
    ax[0].tick_params(
            axis='y', which='both', left=False, right=False,
            labelleft=False)
    ax[0].set_ylim(0, 1.5)
    locs = ax[0].get_xticks()
    time_utc = np.round(np.mod(locs, 1) * 24., 0)
    time_loc = time_utc + utc_offset
    labels = ['{0:.0f}\n{1:.0f}'.format(utc, loc) \
              for (utc, loc) in zip(time_utc, time_loc)]
    ax[0].set_xticklabels(labels)
    ax[0].text(-0.06, -0.03, 'UTC', horizontalalignment='left',
            transform=ax[0].transAxes)
    ax[0].text(-0.06, -0.05, 'UTC{0:+d}'.format(utc_offset),
            horizontalalignment='left', transform=ax[0].transAxes)
    label = block.time_start.value[:10]
    ax[0].set_xlabel(label)

    time_obs = block.time_obs.to_value('h')
    time_slew = block.time_slew.to_value('h')
    time_idle = block.time_idle.to_value('h')
    time_tot = block.duration.to_value('h')
    ax[1].bar(
            [1,2,3], [time_obs,time_slew,time_idle], color=['0.6', '0.7', 'w'],
            linestyle=linestyle, linewidth=linewidth, edgecolor='k')
    offset = 0.01 * ax[1].get_ylim()[1]
    ratio_obs = time_obs / time_tot
    ratio_slew = time_slew / time_tot
    ratio_idle = time_idle / time_tot
    ax[1].text(
            1, time_obs+offset, '{0:0.1f} %'.format(ratio_obs*100.),
            horizontalalignment='center')
    ax[1].text(
            2, time_slew+offset, '{0:0.1f} %'.format(ratio_slew*100.),
            horizontalalignment='center')
    ax[1].text(
            3, time_idle+offset, '{0:0.1f} %'.format(ratio_idle*100.),
            horizontalalignment='center')
    ax[1].set_ylabel('hours')
    labels = ['obs', 'slew', 'idle']
    ax[1].set_xticks([1,2,3])
    ax[1].set_xticklabels(labels)
    plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=45)
    return fig