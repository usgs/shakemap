import os.path
import glob
import logging

import numpy as np

from shakemap.utils.config import get_config_paths
from mapio.gridcontainer import GridHDFContainer
from mapio.grid2d import Grid2D


def get_period_from_imt(imtstr):
    return float(imtstr.replace('SA(', '').replace(')', ''))

def get_generic_amp_factors(sx, myimt):

    # Read the GenericAmpFactors directory in the install directory
    indir, datadir = get_config_paths()
    gaf_dir = os.path.join(indir, "data", "GenericAmpFactors")
    if not os.path.isdir(gaf_dir):
        logging.warn("No GenericAmpFactors directory found.")
        return None
    gaf_files = glob.glob(os.path.join(gaf_dir, '*.hdf'))
    if len(gaf_files) == 0:
        logging.warn("No generic amplification files founs.")
        return None

    gaf = np.zeros_like(sx.lats)

    for gfile in gaf_files:
        thisimt = myimt
        # Get a list of IMTs
        gc = GridHDFContainer.load(gfile)
        contents = gc.getGrids()
        if thisimt == 'PGV' and 'PGV' not in contents:
            logging.warn("Generic Amp Factors: PGV not found in file %s, "
                         "attempting to use SA(1.0)" % (gfile))
            thisimt = 'SA(1.0)'
        if thisimt == 'PGA' and 'PGA' not in contents:
            logging.warn("Generic Amp Factors: PGA not found in file %s, "
                         "attempting to use SA(0.01)" % (gfile))
            thisimt = 'SA(0.01)'

        if thisimt in contents:
            # If imt in IMT list, get the grid
            mygrid, metadata = gc.getGrid(thisimt)
        elif not thisimt.startswith('SA('):
            logging.warn("Generic Amp Factors: IMT %s not found in file %s"
                    % (myimt, gfile))
            mygrid = None
        else:
            # Get the weighted average grid based on the 
            # periods bracketing the input IMT
            mygrid, metadata = _get_average_grid(gc, contents, thisimt)

        if mygrid is None:
            continue

        amp_factors = mygrid.getValue(sx.lats, sx.lons, default=0.0)
        # Sum into output array
        gaf += amp_factors
    return gaf

def _get_average_grid(gc, contents, myimt):
    """
    Given an SA(X) IMT, attempt to find the grids that bracket its
    period and return an interpolated grid that is weighted average
    (weighted by the (log) differeences in period). If the period
    is less than the lowest, or greater than the highest, available
    period, then the closest endpoint grid is returned.

    Args:
        gc (GridHDFContainer): The container holding the amplification
            grids, labeled by IMT string.
        contents (list): A list of the IMTs available in gc.
        myimt (str): The target IMT; must be of type "SA(X)".

    Returns:
        tuple: A grid and its associated metadata.
    """

    #
    # Make a list of the SA IMTs, add the target IMT to the list
    # and then sort by period.
    #
    imt_list = [thisimt for thisimt in contents if thisimt.startswith('SA(')]
    if len(imt_list) == 0:
        logging.warn('Generic Amp Factors: No SA grids in file')
        return None, None
    imt_list.append(myimt)
    imt_list_sorted = sorted(imt_list, key=get_period_from_imt)
    nimt = len(imt_list_sorted)
    ix = imt_list_sorted.index(myimt)
    if ix == 0:
        logging.warn("Generic Amp Factors:IMT %s less than min available "
                     "imt, using %s" % (myimt, imt_list_sorted[1]))
        return gc.getGrid(imt_list_sorted[1])
    elif ix == (nimt - 1):
        logging.warn("Generic Amp Factors:IMT %s greater than max available "
                     "imt, using %s" % (myimt, imt_list_sorted[-2]))
        return gc.getGrid(imt_list_sorted[-2])
    else:
        # Interpolate using (log) period: p1 is the shorter period,
        # p2 is the longer period, and p0 is the target period.
        g1, md1 = gc.getGrid(imt_list_sorted[ix - 1])
        g2, md1 = gc.getGrid(imt_list_sorted[ix + 1])
        p1 = np.log(get_period_from_imt(imt_list_sorted[ix - 1]))
        p2 = np.log(get_period_from_imt(imt_list_sorted[ix + 1]))
        p0 = np.log(get_period_from_imt(myimt))
        w1 = (p2 - p0) / (p2 - p1)
        w2 = 1.0 - w1
        gmean = g1.getData() * w1 + g2.getData() * w2
        return Grid2D(gmean, g1.getGeoDict()), md1

