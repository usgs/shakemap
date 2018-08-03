# stdlib imports
import os.path
import glob
import json
from collections import OrderedDict
import shutil
import logging
import argparse
import inspect

# third party imports
import fiona
from openquake.hazardlib import imt

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths
from shakemap.utils.logging import get_logging_config
from shakelib.plotting.contour import contour
from shakemap.utils.utils import get_object_from_config
from shakelib.utils.containers import ShakeMapOutputContainer
from shakelib.utils.imt_string import oq_to_file

FORMATS = {
    'geojson': ('GeoJSON', 'json')
}

DEFAULT_FILTER_SIZE = 10


class ContourModule(CoreModule):
    """
    contour -- Generate contours of all IMT values from the
                     shake_result.hdf output file.
    """

    command_name = 'contour'
    targets = [r'products/cont_.*\.json']
    dependencies = [('products/shake_result.hdf', True)]

    # supply here a data structure with information about files that
    # can be created by this module.
    contour_page = {'title': 'Ground Motion Contours', 'slug': 'contours'}
    contents = OrderedDict.fromkeys(['miContour', 'pgaContour', 'pgvContour',
                                     'psa[PERIOD]Contour'])
    contents['miContour'] = {'title': 'Intensity Contours',
                             'caption': 'Contours of macroseismic intensity.',
                             'page': contour_page,
                             'formats': [{'filename': 'cont_*mmi.json',
                                          'type': 'application/json'}
                                         ]
                             }
    contents['pgaContour'] = {'title': 'PGA Contours',
                              'caption': 'Contours of [COMPONENT] peak '
                                         'ground acceleration (%g).',
                              'page': contour_page,
                              'formats': [{'filename': 'cont_*pga.json',
                                           'type': 'application/json'}
                                          ]
                              }
    contents['pgvContour'] = {'title': 'PGV Contours',
                              'caption': 'Contours of [COMPONENT] peak '
                                         'ground velocity (cm/s).',
                              'page': contour_page,
                              'formats': [{'filename': 'cont_*pgv.json',
                                           'type': 'application/json'}
                                          ]
                              }
    psacap = 'Contours of [COMPONENT] [FPERIOD] sec 5% damped '\
             'pseudo-spectral acceleration(%g).'
    contents['psa[PERIOD]Contour'] = {
        'title': 'PSA[PERIOD] Contour',
        'page': contour_page,
        'caption': psacap,
        'formats': [{'filename': 'cont_*psa[0-9]p[0-9].json',
                     'type': 'application/json'}
                    ]
    }

    def __init__(self, eventid, filter=None):
        """
        Instantiate a ContourModule class with an event ID.
        """
        self._eventid = eventid
        log_config = get_logging_config()
        log_name = log_config['loggers'].keys()[0]
        self.logger = logging.getLogger(log_name)
        if filter is not None:
            self.filter_size = filter
        else:
            self.filter_size = DEFAULT_FILTER_SIZE

    def execute(self):
        """
        Create contour files for all configured IMT values.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        datafile = os.path.join(datadir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the ShakeMapOutputContainer and extract the data
        container = ShakeMapOutputContainer.load(datafile)

        if container.getDataType() != 'grid':
            raise NotImplementedError('contour module can only contour '
                                      'gridded data, not sets of points')

        # create contour files
        self.logger.debug('Contouring to files...')
        contour_to_files(container, datadir, self.logger, self.filter_size)
        container.close()

    def parseArgs(self, arglist):
        """
        Set up the object to accept the --filter flag.
        """
        parser = argparse.ArgumentParser(
            prog=self.__class__.command_name,
            description=inspect.getdoc(self.__class__))
        parser.add_argument('-f', '--filter', help='Specify the filter '
                            'size (in grid points) for smoothing the '
                            'grids before contouring. Must be a positive'
                            'integer (default=10; use 0 to disable '
                            'filtering).', type=int)
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
        if args.filter is None:
            self.filter_size = DEFAULT_FILTER_SIZE
        else:
            self.filter_size = args.filter
        return args.rem


def contour_to_files(container, output_dir, logger,
                     filter_size=DEFAULT_FILTER_SIZE):
    """
    Generate contours of all IMT values.

    Args:
      container (ShakeMapOutputContainer): ShakeMapOutputContainer with
          ShakeMap output data.
      output_dir (str): Path to directory where output files will be written.
      logger (logging.Logger): Python logging Logger instance.

    Raises:
        LookupError: When configured file format is not supported
    """

    imtlist = container.getIMTs()

    # Right now geojson is all we support; if that changes, we'll have
    # to add a configuration or command-line option
    file_format = 'geojson'
    # open a file for writing
    driver, extension = FORMATS[file_format]
    sa_schema = {
        'geometry': 'MultiLineString',
        'properties': {
            'value': 'float',
            'units': 'str',
            'color': 'str',
            'weight': 'int'
        }
    }
    mmi_schema = {
        'geometry': 'MultiLineString',
        'properties': {
            'value': 'float',
            'units': 'str',
            'color': 'str',
            'weight': 'int'
        }
    }
    crs = {
        'no_defs': True,
        'ellps': 'WGS84',
        'datum': 'WGS84',
        'proj': 'longlat'
    }

    config = container.getConfig()
    gmice = get_object_from_config('gmice', 'modeling', config)
    gmice_imts = gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES
    gmice_pers = gmice.DEFINED_FOR_SA_PERIODS

    for imtype in imtlist:
        fileimt = oq_to_file(imtype)
        oqimt = imt.from_string(imtype)
        component = container.getComponents(imtype)[0]
        if component == 'GREATER_OF_TWO_HORIZONTAL':
            fname = 'cont_%s.%s' % (fileimt, extension)
        else:
            fname = 'cont_%s_%s.%s' % (fileimt, component, extension)
        filename = os.path.join(output_dir, fname)
        if os.path.isfile(filename):
            fpath, fext = os.path.splitext(filename)
            flist = glob.glob(fpath + '.*')
            for fname in flist:
                os.remove(fname)

        if imtype == 'MMI' or not isinstance(oqimt, tuple(gmice_imts)) or \
           (isinstance(oqimt, imt.SA) and oqimt.period not in gmice_pers):
            my_gmice = None
        else:
            my_gmice = gmice

        # fiona spews a warning here when driver is geojson
        # this warning appears to be un-catchable using
        # with warnings.catch_warnings()
        # or
        # logging.captureWarning()
        # or
        # even redirecting stderr/stdout to IO streams
        # not sure where the warning is coming from,
        # but there appears to be no way to stop it...
        with fiona.drivers():
            if imtype == 'MMI':
                selected_schema = mmi_schema
            else:
                selected_schema = sa_schema
            vector_file = fiona.open(
                filename, 'w',
                driver=driver,
                schema=selected_schema,
                crs=crs
            )

            line_strings = contour(container, imtype, component, filter_size,
                                   my_gmice)

            for feature in line_strings:
                vector_file.write(feature)

                # Grab some metadata
            meta = container.getMetadata()
            event_info = meta['input']['event_information']
            mdict = {
                'eventid': event_info['event_id'],
                'longitude': float(event_info['longitude']),
                'latitude': float(event_info['latitude'])
            }

            logger.debug('Writing contour file %s' % filename)
            vector_file.close()

            # Get bounds
            tmp = fiona.open(filename)
            bounds = tmp.bounds
            tmp.close()

            # Read back in to add metadata/bounds
            data = json.load(open(filename))
            data['metadata'] = mdict
            data['bbox'] = bounds
            with open(filename, 'w') as outfile:
                json.dump(data, outfile)

            #####################################
            # Make an extra version of the MMI contour file
            # so that the current web rendering code can find it.
            # Delete this file once everyone has moved to new version
            # of ComCat code.

            if imtype == 'MMI':
                old_file = os.path.join(output_dir, 'cont_mi.json')
                shutil.copy(filename, old_file)
            #####################################
