#stdlib imports
import sys
import os.path
import json
import logging
import glob
import warnings
import io

#third party imports
from shapely.geometry import MultiLineString,mapping
import fiona
import numpy as np
from skimage import measure
from scipy.ndimage.filters import median_filter
from shakelib.utils.containers import OutputContainer
from shakelib.utils.imt_string import oq_to_file, file_to_oq
from configobj import ConfigObj

#local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths,get_logging_config

FORMATS = {'shapefile':('ESRI Shapefile','shp'),
           'geojson':('GeoJSON','json')}

DEFAULT_FILTER_SIZE = 10

class ContourModule(CoreModule):
    """
    contour - Generate contours of all configured IMT values.
    """
    command_name = 'contour'
    def execute(self):
        """Create contour files for all configured IMT values.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not exist.
        """
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        datafile = os.path.join(datadir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the OutputContainer and extract the data
        container = OutputContainer.load(datafile)

        # get the path to the products.conf file, load the config
        config_file = os.path.join(install_path, 'config', 'products.conf')
        config = ConfigObj(config_file)

        # create contour files
        self.logger.info('Contouring to files...')
        contour_to_files(container,config,datadir,self.logger)

def contour(container,imtype,component,intervals=None,
            filter_size=DEFAULT_FILTER_SIZE):
    """Generate contours of a specific IMT and return as a Shapely MultiLineString object.

    Args:
      container (OutputContainer): OutputContainer with ShakeMap output data.
      imtype (str): String containing the name of an Intensity Measure Type 
                    found in container.
      component (str): Intensity Measure component found in container.
      intervals (np.ndarray or None): Array of intervals for IMT, or None.
      filter_size (int): Integer filter (see
                         https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.ndimage.filters.median_filter.html)
    Returns:
        list: List of dictionaries containing two fields:
            - geometry: GeoJSON-like representation of one of the objects in
                    https://toblerity.org/fiona/manual.html#geometry-types
            - properties: Dictionary of properties describing that feature.
    """
    imtdict = container.getIMT(imtype,component)
    gridobj = imtdict['mean']
    grid = gridobj.getData()
    metadata = gridobj.getGeoDict().asDict()
    if imtype == 'MMI':
        sgrid = grid
        fgrid = median_filter(sgrid, size=filter_size)
        units = 'mmi'
    elif imtype == 'PGV':
        sgrid = np.exp(grid)
        fgrid = median_filter(sgrid, size=10)
        units = 'cms'
    else:
        sgrid = np.exp(grid) * 100.0
        fgrid = median_filter(sgrid, size=10)
        units = 'pctg'

    
    if intervals is None:
        interval_type = 'log'
        if imtype == 'MMI':
            interval_type = 'linear'
        intervals = _get_default_intervals(fgrid,interval_type = interval_type)

    lonstart = metadata['xmin']
    latstart = metadata['ymin']
    lonspan = np.abs(metadata['xmax'] - lonstart)
    latspan = np.abs(metadata['ymax'] - latstart)
    nlon = metadata['nx']
    nlat = metadata['ny']

    line_strings = [] #dictionary of MultiLineStrings and props

    for cval in intervals:
        contours = measure.find_contours(fgrid, cval)
        #
        # Convert coords to geographic coordinates; the coordinates
        # are returned in row, column order (i.e., (y, x))
        #
        new_contours = []
        plot_contours = []
        for ic, coords in enumerate(contours): #coords is a line segment
            if len(coords) <= 20: #skipping little contour islands?
                continue

            contours[ic][:, 1] = coords[:, 1] * lonspan / nlon + lonstart
            contours[ic][:, 0] = (nlat - coords[:, 0]) * latspan / nlat + latstart
            plot_contours.append(contours[ic])
            new_contours.append(contours[ic].tolist())

        if len(new_contours):
            mls = MultiLineString(new_contours)
            props = {'value': cval, 'units': units}
            line_strings.append({'geometry':mapping(mls),
                               'properties':props})
    return line_strings

def contour_to_files(container,config,output_dir,logger):
    """Generate contours of all configured IMT values.

    Args:
      container (OutputContainer): OutputContainer with ShakeMap output data.
      config (dict): Product configuration information (from product.conf).
      output_dir (str): Path to directory where output files will be written.
      logger (logging.Logger): Python logging Logger instance.

    Raises:
        LookupError: When configured file format is not supported, or
            when configured IMT is not found in container.
        
    """
    verbose = True
    jsonstr = container.getString('info.json')
    infojson = json.loads(jsonstr)
    event_info = {'event_id': infojson['input']['event_information']['event_id'],
                  'latitude': infojson['input']['event_information']['latitude'],
                  'longitude': infojson['input']['event_information']['longitude']}

    contour_dict = {}
    imtlist = config['products']['contours']['IMTS'].keys()

    file_format =  config['products']['contours']['format']
    #open a file for writing
    if file_format not in FORMATS:
        raise LookupError('File format %s not supported for contours.' % file_format)
    driver,extension = FORMATS[file_format]
    schema = {'geometry':'MultiLineString',
              'properties':{'value':'float',
                            'units':'str'}}
    crs = {'no_defs': True, 'ellps': 'WGS84', 'datum': 'WGS84', 'proj': 'longlat'}
    
    for imtype in imtlist:
        fileimt = oq_to_file(imtype)
        try:
            components = container.getComponents(imtype)
        except LookupError as look_error:
            fmt = 'No IMT called %s in container %s. Skipping.'
            logger.warn(fmt % (imtype,container.getFileName()))
            continue
        imtype_spec = config['products']['contours']['IMTS'][imtype]
        filter_size = int(imtype_spec['filter_size'])
        for component in components:
            fname = '%s_%s.%s' % (fileimt,component,extension)
            filename = os.path.join(output_dir,fname)
            if os.path.isfile(filename):
                fpath,fext = os.path.splitext(filename)
                flist = glob.glob(fpath+'.*')
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
                vector_file = fiona.open(filename,'w',
                                         driver=driver,
                                         schema=schema,
                                         crs=crs)


                intervals = None
                if 'intervals' in imtype_spec:
                    intervals = [float(i) for i in imtype_spec['intervals']]

                line_strings = contour(container,imtype,component,intervals)
                for feature in line_strings:
                  vector_file.write(feature)

                logger.debug('Writing contour file %s' % filename)
                vector_file.close()

def _get_default_intervals(fgrid,interval_type='log'):
    """Get default intervals for any IMT.
    
    Args:
        fgrid (ndarray): Numpy array containing IMT (MMI,PGA, etc.) data.
        interval_type (str): Either 'log' or 'linear'.
    
    Returns:
        ndarray: Numpy array with contour intervals for input IMT.
    """
    if interval_type == 'log':
        #get the range of the input data
        dmin = fgrid.min()
        dmax = fgrid.max()
        #create an array of intervals at 1's and 3's
        intervals = []
        for i in range(-5,5):
            intervals += [10**i,(10**i)*3]
        intervals = np.array(intervals)

        #subtract the intervals from the lowest value in the data
        diffs = dmin - intervals

        #find the index of the largest negative value in diffs
        ibottom = np.where(diffs == diffs[diffs < 0].max())[0][0]
        cmin = intervals[ibottom]

        #find the index
        diffs = dmax - intervals
        diffmin = diffs[diffs > 0].min()
        itop = np.where(diffs == diffmin)[0][0]
        cmax = intervals[itop]
        return intervals[ibottom:itop+1]
    else:
        # Because you dont contour below the smallest value
        gmin = np.ceil(np.min(fgrid))
        # Because you don't contour above the largest value
        gmax = np.floor(np.max(fgrid))
        intervals = np.arange(gmin,gmax+1,1)
        return intervals
