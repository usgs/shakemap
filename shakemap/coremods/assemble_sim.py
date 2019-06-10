# stdlib imports
import os.path
import glob
import re

# third party imports
import numpy as np
import pandas as pd
from configobj import ConfigObj
from mapio.grid2d import Grid2D
from mapio.geodict import GeoDict


# local imports
from .base import CoreModule, Contents
from shakemap.utils.config import get_config_paths

LATLON_COLS = set(['LAT', 'LON'])
XY_COLS = set(['X', 'Y'])

GEO_PROJ_STR = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

IMT_MATCHES = ['MMI',
               'PGA',
               'PGV',
               r'SA\([-+]?[0-9]*\.?[0-9]+\)']


class AssembleSimModule(CoreModule):
    """
    assemble_sim -- Assemble simulation data, configuration information,
    and event.xml into ShakeMap results.
    """

    command_name = 'assemble_sim'
    targets = []
    dependencies = [('event.xml', True),
                    ('simulation.conf', True),
                    ('simulation_*.csv', True)]

    def __init__(self, eventid):
        super(AssembleSimModule, self).__init__(eventid)
        self.contents = Contents(None, None, eventid)

    def execute(self):
        """
        Assemble a simulated ShakeMap from simulated ground motions.

        Raises:
            NotADirectoryError if data directory does not exist.
            FileNotFoundError if simulation.conf or simulation CSV
                could not be found.
            FileExistsError more than one simulation CSV file is found.
        """
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)

        # load the simulation.conf file
        config_file = os.path.join(datadir, 'simulation.conf')
        if not os.path.isfile(config_file):
            raise FileNotFoundError(
                'Could not find simulation config file %s' % config_file)

        # find the spec file for simulation.conf
        config = ConfigObj(infile=config_file)

        # look for the simulation data file
        simfiles = glob.glob(os.path.join(datadir, 'simulation_*.csv'))
        if not len(simfiles):
            raise FileNotFoundError('No simulation data files found.')
        if len(simfiles) > 1:
            raise FileExistsError('Too many simulation data files found.')
        simfile = simfiles[0]

        imtgrids = _get_grids(config, simfile)

        # do stuff with this dictionary of IMT Grid2D objects...


def _get_grids(config, simfile):
    """Create a dictionary of Grid2D objects for each IMT in input CSV file.

    Args:
        config (ConfigObj): Dictionary containing fields:
                            - simulation: (dict)
                              - order = (required) "rows" or "cols"
                              - projection = (optional) Proj4 string
                                defining input X/Y data projection.
                              - nx Number of columns in input grid.
                              - ny Number of rows in input grid.
                              - dx Resolution of columns (if XY, whatever
                                those units are, otherwise decimal degrees).
                              - dy Resolution of rows (if XY, whatever those
                                units are, otherwise decimal degrees).
        simfile (str): Path to a CSV file with columns:
                       - LAT Latitudes for each cell. If irregular, X/Y data
                         will be used.
                       - LON Longitudes for each cell. If irregular, X/Y data
                         will be used.
                       - X Regularized X coordinates for each cell.
                       - Y Regularized Y coordinates for each cell.
                       - H1_<IMT> First horizontal channel for given IMT.
                         Supported IMTs are: PGA, PGV, SA(period).
                       - H2_<IMT> Second horizontal channel for given IMT.
                         Supported IMTs are: PGA, PGV, SA(period).

    Returns:
        dict: Dictionary of IMTs (PGA, PGV, SA(1.0), etc.) and Grid2D
        objects. If XY data was used, these grids are the result of a
        projection/resampling of that XY data back to a regular lat/lon grid.
    """
    row_order = 'C'
    if config['simulation']['order'] != 'rows':
        row_order = 'F'
    dataframe = pd.read_csv(simfile)

    # construct a geodict
    geodict = _get_geodict(dataframe, config)

    # figure out which IMTs we have...
    column_list = []
    for column in dataframe.columns:
        for imtmatch in IMT_MATCHES:
            if re.search(imtmatch, column):
                column_list.append(column)

    # gather up all "channels" for each IMT
    imtdict = {}  # dictionary of imts and a list of columns
    for col in column_list:
        channel, imt = col.split('_')
        if imt in imtdict:
            imtdict[imt].append(col)
        else:
            imtdict[imt] = [col]

    # make a dictionary of Grid2D objects containing max of two "channels"
    # for each IMT
    nrows = geodict.ny
    ncols = geodict.nx
    imtgrids = {}
    for imt, imtcols in imtdict.items():
        icount = len(imtcols)
        if icount != 2:
            raise IndexError(
                'Incorrect number of channels for IMT %s.' % imt)
        channel1, channel2 = imtcols
        maximt = dataframe[[channel1, channel2]].max(axis=1).values
        data = np.reshape(maximt, (nrows, ncols), order=row_order)
        grid = Grid2D(data, geodict)

        # if we need to project data back to geographic, do that here
        if geodict.projection != GEO_PROJ_STR:
            grid = grid.project(GEO_PROJ_STR)

        imtgrids[imt] = grid

    return imtgrids


def _get_geodict(dataframe, config):
    """Get a GeoDict object from input dataframe and simulation config.

    Args:
        dataframe (pandas DataFrame): Input CSV file as a data structure.
        config (ConfigObj): Dictionary with spatial information about CSV file.

    Returns:
        GeoDict: Upper left corner, cell dimensions, rows/cols, and projection.

    """
    nx = config['simulation']['nx']
    ny = config['simulation']['ny']
    dx = config['simulation']['dx']
    dy = config['simulation']['dy']
    row_order = config['simulation']['order'] == 'rows'

    has_latlon = False
    if LATLON_COLS <= set(dataframe.columns):
        has_latlon = True
    has_xy = False
    if XY_COLS <= set(dataframe.columns):
        has_xy = True

    # check for lat/lon regularity
    is_regular = False
    if has_latlon:
        lat = dataframe['LAT'].values
        lon = dataframe['LON'].values
        ymin = lat.min()
        ymax = lat.max()
        xmin = lon.min()
        xmax = lon.max()
        projection = GEO_PROJ_STR
        if row_order:
            lat = np.reshape(lat, (ny, nx), order='C')
            lon = np.reshape(lon, (ny, nx), order='C')
        else:
            lat = np.reshape(lat, (ny, nx), order='F')
            lon = np.reshape(lon, (ny, nx), order='F')
        try:
            np.testing.assert_almost_equal(lat[0][0], lat[0][1])
            is_regular = True
        except AssertionError:
            pass
        if not is_regular and not has_xy:
            raise IndexError(
                'Input simulation files must have *regular* projected X/Y '
                'or lat/lon coordinates.')
    if not is_regular:
        x = dataframe['X'].values
        y = dataframe['X'].values
        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        projection = config['simulation']['projection']

    gd = {'xmin': xmin,
          'xmax': xmax,
          'ymin': ymin,
          'ymax': ymax,
          'nx': nx,
          'ny': ny,
          'dx': dx,
          'dy': dy,
          'projection': projection}
    geodict = GeoDict(gd)

    return geodict

# def unproject_grid2d(grid):
#     geodict = grid.getGeoDict()
#     # establish the input projection
#     input_srs = osr.SpatialReference()
#     input_srs.ImportFromProj4(geodict.projection)

#     output_srs = osr.SpatialReference()
#     output_srs.ImportFromProj4(GEO_PROJ_STR)

#     nrows, ncols = grid._data.shape

#     src_transform = Affine.from_gdal(geodict.xmin -
#                                      geodict.dx / 2.0,
#                                      geodict.dx,
#                                      0.0,  # x rotation, not used by us
#                                      geodict.ymax
#                                      + geodict. dy / 2.0,
#                                      0.0,  # y rotation, not used by us
#                                      # their dy is negative
#                                      -1 * geodict.dy)
