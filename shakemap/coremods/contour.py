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
from shapely.geometry import MultiLineString
from shapely.geometry import mapping
import fiona
import numpy as np
from skimage import measure
from scipy.ndimage.filters import median_filter
from shakelib.utils.containers import ShakeMapOutputContainer
from shakelib.utils.imt_string import oq_to_file

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths, get_logging_config
from impactutils.colors.cpalette import ColorPalette

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

    # supply here a data structure with information about files that
    # can be created by this module.
    contour_page = {'title': 'Ground Motion Contours', 'slug': 'contours'}
    contents = OrderedDict.fromkeys(['miContour', 'pgaContour', 'pgvContour',
                                     'psa[PERIOD]Contour'])
    contents['miContour'] = {'title': 'Intensity Contours',
                             'caption': 'Contours of macroseismic intensity.',
                             'page': contour_page,
                             'formats': [{'filename': 'cont_*MMI.json',
                                         'type': 'application/json'}
                                        ]
    }
    contents['pgaContour'] = {'title': 'PGA Contours',
                              'caption': 'Contours of [COMPONENT] peak '\
                                         'ground acceleration (%g).',
                              'page': contour_page,
                              'formats': [{'filename': 'cont_*PGA.json',
                                          'type': 'application/json'}
                                         ]
    }
    contents['pgvContour'] = {'title': 'PGV Contours',
                              'caption': 'Contours of [COMPONENT] peak '\
                                         'ground velocity (cm/s).',
                              'page': contour_page,
                              'formats': [{'filename': 'cont_*PGV.json',
                                          'type': 'application/json'}
                                         ]
    }
    psacap = 'Contours of [COMPONENT] [FPERIOD] sec 5% damped '\
             'pseudo-spectral acceleration(%g).'
    contents['psa[PERIOD]Contour'] = {
            'title': 'PSA[PERIOD] Contour',
            'page': contour_page,
            'caption': psacap,
            'formats': [{'filename':'cont_*PSA[0-9]p[0-9].json',
                         'type':'application/json'}
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

    def parseArgs(self, arglist):
        """
        Set up the object to accept the --filter flag.
        """
        parser = argparse.ArgumentParser(prog=self.__class__.command_name,
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


def contour(container, imtype, component, filter_size):
    """
    Generate contours of a specific IMT and return as a Shapely
    MultiLineString object.

    Args:
        container (ShakeMapOutputContainer): ShakeMapOutputContainer
            with ShakeMap output data.
        imtype (str): String containing the name of an Intensity
            Measure Type found in container.
        component (str): Intensity Measure component found in container.
        filter_size (int): Integer filter (see
            https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.ndimage.filters.median_filter.html)
    Returns:
        list: List of dictionaries containing two fields

                - geometry: GeoJSON-like representation of one of the objects
                  in https://toblerity.org/fiona/manual.html#geometry-types
                - properties: Dictionary of properties describing that
                  feature.

    Raises:
        NotImplementedError -- if the user attempts to contour a data file
            with sets of points rather than grids.
    """  # noqa
    intensity_colormap = ColorPalette.fromPreset('mmi')
    imtdict = container.getIMTGrids(imtype, component)
    gridobj = imtdict['mean']
    grid = gridobj.getData()
    metadata = gridobj.getGeoDict().asDict()
    if imtype == 'MMI':
        sgrid = grid
        units = 'mmi'
    elif imtype == 'PGV':
        sgrid = np.exp(grid)
        units = 'cms'
    else:
        sgrid = np.exp(grid) * 100.0
        units = 'pctg'
    if filter_size > 0:
        fgrid = median_filter(sgrid, size=filter_size)
    else:
        fgrid = sgrid

    interval_type = 'log'
    if imtype == 'MMI':
        interval_type = 'linear'
    intervals = getContourLevels(
        np.min(fgrid), np.max(fgrid), itype=interval_type)

    lonstart = metadata['xmin']
    latstart = metadata['ymin']
    lonspan = np.abs(metadata['xmax'] - lonstart)
    latspan = np.abs(metadata['ymax'] - latstart)
    nlon = metadata['nx']
    nlat = metadata['ny']

    line_strings = []  # dictionary of MultiLineStrings and props

    for cval in intervals:
        contours = measure.find_contours(fgrid, cval)
        #
        # Convert coords to geographic coordinates; the coordinates
        # are returned in row, column order (i.e., (y, x))
        #
        new_contours = []
        plot_contours = []
        for ic, coords in enumerate(contours):  # coords is a line segment
            if len(coords) <= 20:  # skipping little contour islands?
                continue

            mylons = coords[:, 1] * lonspan / nlon + lonstart
            mylats = (nlat - coords[:, 0]) * latspan / nlat + latstart
            contours[ic][:, 0] = mylons[:]
            contours[ic][:, 1] = mylats[:]
            plot_contours.append(contours[ic])
            new_contours.append(contours[ic].tolist())

        if len(new_contours):
            mls = MultiLineString(new_contours)
            props = {
                'value': cval,
                'units': units
            }
            if imtype == 'MMI':
                color_array = np.array(intensity_colormap.getDataColor(cval))
                color_rgb = np.array(
                    color_array[0:3] * 255, dtype=int).tolist()
                props['color'] = '#%02x%02x%02x' % tuple(color_rgb)
                if (cval * 2) % 2 == 1:
                    props['weight'] = 4
                else:
                    props['weight'] = 2
            line_strings.append(
                {
                    'geometry': mapping(mls),
                    'properties': props
                }
            )
    return line_strings


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
            'units': 'str'
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

    for imtype in imtlist:
        fileimt = oq_to_file(imtype)
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

            line_strings = contour(
                container,
                imtype,
                component,
                filter_size
            )

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

def getContourLevels(dmin, dmax, itype='log'):
    """
    Get contour levels given min/max values and desired IMT.

    Use itype='log' for any IMT that is logarithmically distributed, such as
    PGA, PGV, and Sa. Linear for MMI.

    Args:
        dmin (float): Minimum value of data to contour.
        dmax (float): Maximum value of data to contour.
        itype (str): Interval type; default is 'log', anythign else
            indicates linear intervals.

    Returns:
        ndarray: Numpy array of contour levels.

    """
    if itype == 'log':
        # Within-decade label values
        dec_inc = np.array([1, 2, 5], dtype=float)

        # Upper and lower decades
        lower_dec = np.floor(np.log10(dmin))
        upper_dec = np.ceil(np.log10(dmax))

        # Array of decades
        decades = np.arange(lower_dec, upper_dec + 1)

        # Construct levels
        levels = np.concatenate([np.power(10, d) * dec_inc for d in decades])
        levels = levels[(levels < dmax) & (levels > dmin)]
    else:
        # MMI contours are every 0.5 units
        levels = np.arange(
            np.ceil(dmin * 2) / 2,
            np.floor(dmax * 2) / 2 + 0.5,
            0.5
        )
    return levels

