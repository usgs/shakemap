# stdlib imports
import os.path
import sys
import glob
import re
import shutil
import datetime
import argparse
import inspect

# third party imports
import numpy as np
import pandas as pd
from configobj import ConfigObj
from mapio.grid2d import Grid2D
import matplotlib.pyplot as plt
from mapio.geodict import GeoDict

# local imports
from .base import CoreModule, Contents
from shakemap.utils.config import (get_config_paths,
                                   get_model_config)
from shakelib.utils.containers import ShakeMapInputContainer
from impactutils.io.smcontainers import ShakeMapOutputContainer
from shakemap.utils.amps import AmplitudeHandler
from shakelib.rupture import constants

LATLON_COLS = set(['LAT', 'LON'])
XY_COLS = set(['X', 'Y'])

GEO_PROJ_STR = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

IMT_MATCHES = ['MMI',
               'PGA',
               'PGV',
               r'SA\([-+]?[0-9]*\.?[0-9]+\)']

SAVE_FILE = '.saved'


class ModelSimModule(CoreModule):
    """
    model_sim -- Assemble simulation data, configuration information,
    and event.xml into ShakeMap results.
    """

    command_name = 'model_sim'
    targets = []
    dependencies = [('event.xml', True),
                    ('simulation.conf', True),
                    ('simulation_*.csv', True)]

    apply_site_amp = False

    rock_vs30 = 760.0

    def __init__(self, eventid, comment=None):
        super(ModelSimModule, self).__init__(eventid)
        self.contents = Contents(None, None, eventid)
        if comment is not None:
            self.comment = comment

    def parseArgs(self, arglist):
        """
        Allow for options:
            -s, --apply_site_amp: Apply site amplification.
        """
        parser = argparse.ArgumentParser(
            prog=self.__class__.command_name,
            description=inspect.getdoc(self.__class__))
        parser.add_argument('-s', '--apply_site_amp', action='store_true',
                            help='Apply site amplification based upon '
                                 'the configured GMPE set.')
        #
        # This line should be in any modules that overrides this
        # one. It will collect up everything after the current
        # modules options in args.rem, which should be returned
        # by this function. Note: doing parser.parse_known_args()
        # will not work as it will suck up any later modules'
        # options that are the same as this one's.
        #
        parser.add_argument('rem', nargs=argparse.REMAINDER,
                            help=argparse.SUPPRESS)
        args = parser.parse_args(arglist)
        if args.apply_site_amp:
            self.apply_site_amp = True
        return args.rem

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

        eventxml = os.path.join(datadir, 'event.xml')
        self.logger.debug('Looking for event.xml file...')
        if not os.path.isfile(eventxml):
            raise FileNotFoundError('%s does not exist.' % eventxml)

        # Prompt for a comment string if none is provided on the command line
        if self.comment is None:
            if sys.stdout is not None and sys.stdout.isatty():
                self.comment = input(
                    'Please enter a comment for this version.\n'
                    'comment: ')
            else:
                self.comment = ''

        # load the simulation.conf file
        config_file = os.path.join(datadir, 'simulation.conf')
        if not os.path.isfile(config_file):
            raise FileNotFoundError(
                'Could not find simulation config file %s' % config_file)

        # find the spec file for simulation.conf
        config = ConfigObj(infile=config_file)

        model_global_config = get_model_config(install_path, datadir,
                                               self.logger)
        model_config = model_global_config.dict()

        self.logger.debug('Looking for rupture files...')
        # look for geojson versions of rupture files
        rupturefile = os.path.join(datadir, 'rupture.json')
        if not os.path.isfile(rupturefile):
            # failing any of those, look for text file versions
            rupturefiles = glob.glob(os.path.join(datadir, '*_fault.txt'))
            rupturefile = None
            if len(rupturefiles):
                rupturefile = rupturefiles[0]

        # find any source.txt or moment.xml files
        momentfile = os.path.join(datadir, 'moment.xml')
        sourcefile = os.path.join(datadir, 'source.txt')
        if not os.path.isfile(sourcefile):
            sourcefile = None
        if not os.path.isfile(momentfile):
            momentfile = None

        #
        # Clear away results from previous runs
        #
        products_path = os.path.join(datadir, 'products')
        if os.path.isdir(products_path):
            shutil.rmtree(products_path, ignore_errors=True)
        pdl_path = os.path.join(datadir, 'pdl')
        if os.path.isdir(pdl_path):
            shutil.rmtree(pdl_path, ignore_errors=True)

        # Look for any .transferred file and delete it
        save_file = os.path.join(datadir, SAVE_FILE)
        if os.path.isfile(save_file):
            os.remove(save_file)

        #
        # Sort out the version history. Get the most recent backup file and
        # extract the existing history. Then add a new line for this run.
        #
        timestamp = datetime.datetime.utcnow().strftime('%FT%TZ')
        originator = config['system']['source_network']
        backup_dirs = sorted(
            glob.glob(os.path.join(datadir, '..', 'backup*')),
            reverse=True)
        if len(backup_dirs):
            #
            # Backup files exist so find the latest one and extract its
            # history, then add a new line that increments the version
            #
            bu_file = os.path.join(backup_dirs[0], 'shake_data.hdf')
            bu_ic = ShakeMapInputContainer.load(bu_file)
            history = bu_ic.getVersionHistory()
            bu_ic.close()
            version = int(
                backup_dirs[0].replace(
                    os.path.join(datadir, '..', 'backup'), ''))
            version += 1
            new_line = [timestamp, originator, version, self.comment]
            history['history'].append(new_line)
        elif os.path.isfile(os.path.join(datadir, 'shake_data.hdf')):
            #
            # No backups are available, but there is an existing shake_data
            # file. Extract its history and update the timestamp and
            # source network (but leave the version alone).
            # If there is no history, just start a new one with version 1
            #
            bu_file = os.path.join(datadir, 'shake_data.hdf')
            bu_ic = ShakeMapInputContainer.load(bu_file)
            history = bu_ic.getVersionHistory()
            bu_ic.close()
            if 'history' in history:
                new_line = [timestamp, originator, history['history'][-1][2],
                            self.comment]
                history['history'][-1] = new_line
            else:
                history = {'history': []}
                new_line = [timestamp, originator, 1, self.comment]
                history['history'].append(new_line)
        else:
            #
            # No backup and no existing file. Make this version 1
            #
            history = {'history': []}
            new_line = [timestamp, originator, 1, self.comment]
            history['history'].append(new_line)

        hdf_file = os.path.join(datadir, 'shake_data.hdf')

        self.logger.debug('Creating input container...')
        shake_data = ShakeMapInputContainer.createFromInput(
            hdf_file,
            config,
            eventxml,
            history,
            rupturefile=rupturefile,
            sourcefile=sourcefile,
            momentfile=momentfile,
            datafiles=[])
        self.logger.debug('Created HDF5 input container in %s' %
                          shake_data.getFileName())
        ah = AmplitudeHandler(install_path, data_path)
        event = ah.getEvent(self._eventid)
        if event is None:
            origin = shake_data.getRuptureObject().getOrigin()
            event = {'id': self._eventid,
                     'netid': origin.netid,
                     'network': origin.network,
                     'time': origin.time.strftime(constants.TIMEFMT),
                     'lat': origin.lat,
                     'lon': origin.lon,
                     'depth': origin.depth,
                     'mag': origin.mag,
                     'locstring': origin.locstring}
            ah.insertEvent(event)
        shake_data.close()

        self.gmice = get_object_from_config('gmice', 'modeling', model_config)

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
    geodict, top_down = _get_geodict(dataframe, config)

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
        if not top_down:
            data = np.flipud(data)
        grid = Grid2D(data, geodict)

        # if we need to project data back to geographic, do that here
        if geodict.projection != GEO_PROJ_STR:
            grid2 = grid.project(GEO_PROJ_STR)

        # remove any nan's, maximizing the resulting area of good data
        grid3 = _trim_grid(grid2)

        imtgrids[imt] = grid3

    return imtgrids


def _trim_grid(ingrid):
    outgrid = Grid2D.copyFromGrid(ingrid)
    while np.isnan(outgrid._data).any():
        nrows, ncols = outgrid._data.shape
        top = outgrid._data[0, :]
        bottom = outgrid._data[-1, :]
        left = outgrid._data[:, 0]
        right = outgrid._data[:, -1]
        ftop = np.isnan(top).sum() / ncols
        fbottom = np.isnan(bottom).sum() / ncols
        fleft = np.isnan(left).sum() / nrows
        fright = np.isnan(right).sum() / nrows
        side = np.argmax([ftop, fbottom, fleft, fright])
        gdict = outgrid.getGeoDict().asDict()
        if side == 0:  # removing top row
            outgrid._data = outgrid._data[1:, :]
            gdict['ymax'] -= gdict['dy']
            gdict['ny'] -= 1
        elif side == 1:  # removing bottom row
            outgrid._data = outgrid._data[0:-1, :]
            gdict['ymin'] += gdict['dy']
            gdict['ny'] -= 1
        elif side == 2:  # removing left column
            outgrid._data = outgrid._data[:, 1:]
            gdict['xmin'] += gdict['dx']
            gdict['nx'] -= 1
        elif side == 3:  # removing right column
            outgrid._data = outgrid._data[:, 0:-1]
            gdict['xmax'] -= gdict['dx']
            gdict['nx'] -= 1
        geodict = GeoDict(gdict)
        outgrid = Grid2D(data=outgrid._data, geodict=geodict)

    return outgrid


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

    # this lets the user know whether the data starts at the top or the bottom
    top_down = True
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
        if lat[1][0] > lat[0][0]:
            top_down = False
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
        y = dataframe['Y'].values
        y2 = np.reshape(y, (ny, nx), order='C')
        if y2[1][0] > y2[0][0]:
            top_down = False
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

    return (geodict, top_down)

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
