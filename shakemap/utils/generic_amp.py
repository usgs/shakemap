import os.path
import glob
import logging

import numpy as np

from shakemap.utils.config import get_config_paths
from mapio.gridcontainer import GridHDFContainer
from mapio.grid2d import Grid2D


def get_period_from_imt(imtstr):
    return float(imtstr.replace("SA(", "").replace(")", ""))


def get_generic_amp_factors(sx, myimt):
    """
    Returns an array of generic amplification factors the same shape
    as the lons and lats members of sx. The amplification factors are
    pulled from the grids in the operator's configured directory:
    '<install>/data/GenericAmpFactors'. Any HDF files in that directory
    will be considerd amplification files and will be examined for
    factors that apply to the coordinates in sx. The factors are
    assumed to be in natural log units and will be added to the output
    of the GMPE (or IPE). For coordinates outside the available grids,
    zero will be returned. If more than one file contains amplifications
    for any of the coordinates (i.e., overlapping grids), then the
    factors will be summed together. It is the operator's responsibility
    to ensure the proper behavior is observed for the selected output
    IMTs considering:

        - If 'myimt' is 'PGA' and there is no PGA grid in the HDF
          file, the 'myimt' will be set to 'SA(0.01)' for that file.
        - If 'myimt' is 'PGV' and there is no PGV grid in the HDF
          file, the 'myimt' will be set to 'SA(1.0)' for that file.
        - If 'myimt' is spectral acceleration (i.e., 'SA(x)' where
          'x' is the period), and SA of that period is not found in the
          HDF file, the function will first attempt to interpolate the
          grids of the next lower and next greater periods found in
          the file. The interpolation is done as a weighted average
          of the grids with the weights being defined assigned by
          the log difference in periods. If the period of 'myimt' is
          less that the shortest period in the file, the grid for
          the shortest period is used. If the period of 'myimt' is
          greater that the longest period in the file, the grid for
          the longest period is used.
        - Interpolation in geographic space is nearest neighbor.
        - Coordinates that fall outside the grid bounds of any
          given file are assigned an amplification factor of zero.

    Args:
        sx (Sites Context): An OpenQuake sites context specifying the
            coordinates of interest.
        myimt (str): A string representing an OpenQuake IMT.

    Returns:
        array: An array of generic amplification factors corresponding
        to the coordinates specified by sx.
    """

    # Read the GenericAmpFactors directory in the install directory
    indir, _ = get_config_paths()
    gaf_dir = os.path.join(indir, "data", "GenericAmpFactors")
    if not os.path.isdir(gaf_dir):
        logging.warning("No GenericAmpFactors directory found.")
        return None
    gaf_files = glob.glob(os.path.join(gaf_dir, "*.hdf"))
    if len(gaf_files) == 0:
        logging.warning("No generic amplification files found.")
        return None

    gaf = np.zeros_like(sx.lats)

    for gfile in gaf_files:
        thisimt = myimt
        # Get a list of IMTs
        gc = GridHDFContainer.load(gfile)
        contents = gc.getGrids()
        if thisimt == "PGV" and "PGV" not in contents:
            logging.warning(
                "Generic Amp Factors: PGV not found in file %s, "
                "attempting to use SA(1.0)" % (gfile)
            )
            thisimt = "SA(1.0)"
        if thisimt == "PGA" and "PGA" not in contents:
            logging.warning(
                "Generic Amp Factors: PGA not found in file %s, "
                "attempting to use SA(0.01)" % (gfile)
            )
            thisimt = "SA(0.01)"

        if thisimt in contents:
            # If imt in IMT list, get the grid
            mygrid, _ = gc.getGrid(thisimt)
        elif not thisimt.startswith("SA("):
            logging.warning(
                f"Generic Amp Factors: IMT {myimt} not found in file {gfile}"
            )
            mygrid = None
        else:
            # Get the weighted average grid based on the
            # periods bracketing the input IMT
            mygrid, metadata = _get_average_grid(gc, contents, thisimt)

        gc.close()
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
    imt_list = [thisimt for thisimt in contents if thisimt.startswith("SA(")]
    if len(imt_list) == 0:
        logging.warning("Generic Amp Factors: No SA grids in file")
        return None, None
    imt_list.append(myimt)
    imt_list_sorted = sorted(imt_list, key=get_period_from_imt)
    nimt = len(imt_list_sorted)
    ix = imt_list_sorted.index(myimt)
    if ix == 0:
        logging.warning(
            "Generic Amp Factors:IMT %s less than min available "
            "imt, using %s" % (myimt, imt_list_sorted[1])
        )
        return gc.getGrid(imt_list_sorted[1])
    elif ix == (nimt - 1):
        logging.warning(
            "Generic Amp Factors:IMT %s greater than max "
            "available imt, using %s" % (myimt, imt_list_sorted[-2])
        )
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
