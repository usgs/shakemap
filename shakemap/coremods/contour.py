# stdlib imports
import os.path
import glob
import json

# third party imports
from shapely.geometry import MultiLineString
from shapely.geometry import mapping
import fiona
import numpy as np
from skimage import measure
from scipy.ndimage.filters import median_filter
from shakelib.utils.containers import ShakeMapOutputContainer
from shakelib.utils.imt_string import oq_to_file
from configobj import ConfigObj

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths
from impactutils.colors.cpalette import ColorPalette

FORMATS = {
    'shapefile': ('ESRI Shapefile', 'shp'),
    'geojson': ('GeoJSON', 'json')
}

DEFAULT_FILTER_SIZE = 10


class ContourModule(CoreModule):
    """
    contour -- Generate contours of all configured IMT values from the
                     shake_result.hdf output file.
    """

    command_name = 'contour'

    def execute(self):
        """Create contour files for all configured IMT values.

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

        # get the path to the products.conf file, load the config
        config_file = os.path.join(install_path, 'config', 'products.conf')
        config = ConfigObj(config_file)

        # create contour files
        self.logger.info('Contouring to files...')
        contour_to_files(container, config, datadir, self.logger)


def contour(container, imtype, component,
            filter_size=DEFAULT_FILTER_SIZE):
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
    """
    intensity_colormap = ColorPalette.fromPreset('mmi')
    if container.getDataType() != 'grid':
        raise NotImplementedError('contour module can only contour '
                                  'gridded data, not sets of points')
    imtdict = container.getIMTGrids(imtype, component)
    gridobj = imtdict['mean']
    grid = gridobj.getData()
    metadata = gridobj.getGeoDict().asDict()
    if imtype == 'MMI':
        sgrid = grid
        fgrid = median_filter(sgrid, size=filter_size)
        units = 'mmi'
    elif imtype == 'PGV':
        sgrid = np.exp(grid)
        fgrid = median_filter(sgrid, size=filter_size)
        units = 'cms'
    else:
        sgrid = np.exp(grid) * 100.0
        fgrid = median_filter(sgrid, size=filter_size)
        units = 'pctg'

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


def contour_to_files(container, config, output_dir, logger):
    """
    Generate contours of all configured IMT values.

    Args:
      container (ShakeMapOutputContainer): ShakeMapOutputContainer with
          ShakeMap output data.
      config (dict): Product configuration information (from product.conf).
      output_dir (str): Path to directory where output files will be written.
      logger (logging.Logger): Python logging Logger instance.

    Raises:
        LookupError: When configured file format is not supported, or
            when configured IMT is not found in container.

    """

    imtlist = config['products']['contours']['IMTS'].keys()

    file_format = config['products']['contours']['format']
    # open a file for writing
    if file_format not in FORMATS:
        raise LookupError(
            'File format %s not supported for contours.' % file_format)
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
        try:
            components = container.getComponents(imtype)
        except LookupError as look_error:
            fmt = 'No IMT called %s in container %s. Skipping.'
            logger.warn(fmt % (imtype, container.getFileName()))
            continue
        imtype_spec = config['products']['contours']['IMTS'][imtype]
        filter_size = int(imtype_spec['filter_size'])
        for component in components:
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


def getContourLevels(dmin, dmax, itype='log'):
    """Get contour levels given min/max values and desired IMT.

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
