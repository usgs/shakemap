"""
Collect configuration, station data, finite fault data, etc., into
an InputContainer and write it out as shake_data.hdf.
"""

# stdlib imports
import argparse
import inspect
import os.path
import glob
import datetime
import shutil
import sys
import re

# third party imports
from configobj import ConfigObj
from validate import Validator
import numpy as np
import pandas as pd
from mapio.grid2d import Grid2D
from mapio.geodict import GeoDict

# local imports
from .base import CoreModule
from shakelib.utils.containers import ShakeMapInputContainer
from shakemap.utils.config import (get_config_paths,
                                   get_configspec,
                                   config_error,
                                   get_model_config,
                                   path_macro_sub)
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


class AssembleModule(CoreModule):
    """
    assemble -- Assemble ShakeMap input data into the shake_data.hdf input
                      file.
    """

    command_name = 'assemble'
    targets = [r'shake_data\.hdf']
    dependencies = [('event.xml', True), ('*_dat.xml', False),
                    ('*_fault.txt', False), ('rupture.json', False),
                    ('source.txt', False), ('model.conf', False),
                    ('model_select.conf', False)]
    configs = ['gmpe_sets.conf', 'model.conf', 'modules.conf']

    def __init__(self, eventid, comment=None):
        """
        Instantiate a CoreModule class with an event ID.
        """
        super(AssembleModule, self).__init__(eventid)
        if comment is not None:
            self.comment = comment

    def execute(self):
        """
        Assemble ShakeMap input data and write and ShakeMapInputContainer named
        shake_data.hdf in the event's 'current' directory.

        Raises:
            NotADirectoryError: When the event data directory does not
                exist.
            FileNotFoundError: When the the event's event.xml file does
                not exist.
            RuntimeError: When there are problems parsing the configuration.
            ValidateError: When there are configuration items missing or mis-
                configured.
        """

        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
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
        # Get the combined model config file
        #
        global_config = get_model_config(install_path, datadir, self.logger)

        global_data_path = os.path.join(os.path.expanduser('~'),
                                        'shakemap_data')
        #
        # If there is a prediction_location->file file, then we need
        # to expand any macros; this could have the event ID, so we
        # can't just use the file_type handler in the configspec
        #
        if 'file' in global_config['interp']['prediction_location']:
            loc_file = global_config['interp']['prediction_location']['file']
            if loc_file and loc_file != 'None':      # 'None' is a string here
                loc_file = path_macro_sub(loc_file, ip=install_path,
                                          dp=data_path, gp=global_data_path,
                                          ei=self._eventid)
                if not os.path.isfile(loc_file):
                    raise FileNotFoundError("prediction file '%s' is not "
                                            "a valid file" % loc_file)
                global_config['interp']['prediction_location']['file'] = \
                    loc_file

        config = global_config.dict()

        self.logger.debug('Looking for data files...')
        datafiles = glob.glob(os.path.join(datadir, '*_dat.xml'))
        if os.path.isfile(os.path.join(datadir, 'stationlist.xml')):
            datafiles.append(os.path.join(datadir, 'stationlist.xml'))
        datafiles += glob.glob(os.path.join(datadir, '*_dat.json'))
        if os.path.isfile(os.path.join(datadir, 'stationlist.json')):
            datafiles.append(os.path.join(datadir, 'stationlist.json'))

        self.logger.debug('Looking for rupture files...')
        # look for geojson versions of rupture files
        rupturefile = os.path.join(datadir, 'rupture.json')
        if not os.path.isfile(rupturefile):
            # failing any of those, look for text file versions
            rupturefiles = glob.glob(os.path.join(datadir, '*_fault.txt'))
            rupturefile = None
            if len(rupturefiles):
                rupturefile = rupturefiles[0]

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
            datafiles=datafiles)
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

        # Look for grids of simulated data
        simfiles = glob.glob(os.path.join(datadir, 'simulation_*.csv'))
        if len(simfiles) > 1:
            raise FileExistsError('Too many simulation data files found.')
        elif len(simfiles):
            # load the simulation.conf file
            config_file = os.path.join(datadir, 'simulation.conf')
            if not os.path.isfile(config_file):
                raise FileNotFoundError(
                    'Could not find simulation config file %s' % config_file)

            # find the spec file for simulation.conf
            sim_config_spec = get_configspec('simulation')
            sim_config = ConfigObj(infile=config_file,
                                   configspec=sim_config_spec)
            results = sim_config.validate(Validator())
            if not isinstance(results, bool) or not results:
                config_error(global_config, results)

            simfile = simfiles[0]
            imtgrids = _get_grids(sim_config, simfile)
            for imtstr in imtgrids:
                metadata = imtgrids[imtstr].getGeoDict().asDict()
                datagrid = imtgrids[imtstr].getData()
                shake_data.setArray(['simulations'], imtstr, datagrid,
                                    metadata=metadata)

        shake_data.close()

    def parseArgs(self, arglist):
        """
        Set up the object to accept the --comment flag.
        """
        parser = argparse.ArgumentParser(
            prog=self.__class__.command_name,
            description=inspect.getdoc(self.__class__))
        parser.add_argument('-c', '--comment', help='Provide a comment for '
                            'this version of the ShakeMap. If the comment '
                            'has spaces, the string should be quoted (e.g., '
                            '--comment "This is a comment.")')
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
        self.comment = args.comment
        return args.rem


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
        data = np.log(np.reshape(maximt, (nrows, ncols), order=row_order))
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
