# stdlib imports
import os.path
import zipfile
import tempfile
from shutil import copyfile
import concurrent.futures as cf
from collections import OrderedDict
from functools import partial

# third party imports
import fiona
from impactutils.io.smcontainers import ShakeMapOutputContainer
import numpy as np
from configobj import ConfigObj
from openquake.hazardlib import imt as OQIMT

# local imports
from .base import CoreModule, Contents
from shakemap.utils.config import (get_data_path,
                                   get_config_paths,
                                   get_configspec,
                                   get_custom_validator,
                                   config_error)
from shakemap.utils.utils import get_object_from_config
from shakelib.plotting.contour import contour
from shakemap.c.pcontour import pcontour
from shakelib.utils.imt_string import oq_to_file


class ShapeModule(CoreModule):
    """
    shape -- Generate shape files of the ground motion parameters
    """

    command_name = 'shape'
    targets = [r'products/shape\.zip']
    dependencies = [('products/shake_result.hdf', True)]

    def __init__(self, eventid):
        """
        Instantiate a ShapeModule class with an event ID.
        """
        super(ShapeModule, self).__init__(eventid)
        self.contents = Contents('GIS Shape Files', 'shape files', eventid)

    def execute(self):
        """
        Create shape files.

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
            raise NotImplementedError('shape module can only contour '
                                      'gridded data, not sets of points')

        config_file = os.path.join(install_path, 'config', 'products.conf')
        spec_file = get_configspec('products')
        validator = get_custom_validator()
        config = ConfigObj(config_file, configspec=spec_file)
        results = config.validate(validator)
        if not isinstance(results, bool) or not results:
            config_error(config, results)

        max_workers = config['products']['mapping']['max_workers']
        method = config['products']['shape']['method']

        create_polygons(container, datadir, self.logger, max_workers,
                        method=method)

        container.close()

        self.contents.addFile('shakemap_shapefiles', 'ShakeMap Shape Files',
                              'Shape Files.',
                              'shape.zip', 'application/zip')


def create_polygons(container, datadir, logger, max_workers, method='pcontour'):
    """ Generates a set of closed polygons (with or without holes) using the
    specified method (either pcontour or skimage), and uses fiona to convert
    the resulting GeoJSON objects into ESRI-style shape files which are then
    zipped into an archive along with .prj, .lyr, and metadata .xml files. A
    warning will be emitted if .lyr, or .xml files cannot be found for the
    ground motion parameter in question.

    Args:
        container (ShakeMapOutputContainer): An open ShakeMap output
            container object.
        datadir (str): The products directory for the event in question.
        logger (logger): This module's logger object.
        method (str): Contouring implementation to use (either 'pcontour' or
            'skimage')

    Returns:
        (nothing): Nothing.
    """

    # gmice info for shakelib.plotting.contour
    config = container.getConfig()
    gmice = get_object_from_config('gmice', 'modeling', config)
    gmice_imts = gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES
    gmice_pers = gmice.DEFINED_FOR_SA_PERIODS

    component = list(container.getComponents())[0]
    imts = container.getIMTs(component)

    if method == 'pcontour':
        schema = {'properties': OrderedDict([('AREA', 'float:13.3'),
                                             ('PERIMETER', 'float:14.3'),
                                             ('PGAPOL_', 'int:12'),
                                             ('PGAPOL_ID', 'int:12'),
                                             ('GRID_CODE', 'int:12'),
                                             ('PARAMVALUE', 'float:14.4')]),
                  'geometry': 'Polygon'}
    elif method == 'skimage':
        schema = {'properties': OrderedDict([('value', 'float:2.1'),
                                             ('units', 'str'),
                                             ('color', 'str'),
                                             ('weight', 'float:13.3')]),
                  'geometry': 'MultiLineString'}
    else:
        raise ValueError('Unknown contouring method {}'.format(method))

    smdata = os.path.join(get_data_path(), 'gis')
    # Make a directory for the files to live in prior to being zipped
    alist = []
    with tempfile.TemporaryDirectory(dir=datadir) as tdir:
        for imt in imts:
            gdict = container.getIMTGrids(imt, component)
            fgrid = gdict['mean']
            if imt == 'MMI':
                fname = 'mi'
            elif imt == 'PGV':
                fname = 'pgv'
            else:
                fname = oq_to_file(imt)

            if method == 'pcontour':
                my_gmice = None
                if imt == 'MMI':
                    contour_levels = np.arange(0.1, 10.2, 0.2)
                elif imt == 'PGV':
                    fgrid = np.exp(fgrid)
                    cont_max = np.ceil(np.max(fgrid)) + 2.0
                    contour_levels = np.arange(1.0, cont_max, 2.0)
                    if contour_levels.size == 0:
                        contour_levels = np.array([1.0])
                else:
                    fgrid = np.exp(fgrid)
                    cont_max = (np.ceil(100 * np.max(fgrid)) + 2.0) / 100.0
                    contour_levels = np.arange(0.01, cont_max, 0.02)
                    if contour_levels.size == 0:
                        contour_levels = np.array([0.01])
            else:
                # skimage method chooses its own levels
                contour_levels = None
                # but wants gmice info
                oqimt = OQIMT.from_string(imt)
                if imt == 'MMI' or not isinstance(oqimt, tuple(gmice_imts)) or \
                   (isinstance(oqimt, OQIMT.SA) and oqimt.period not in gmice_pers):
                    my_gmice = None
                else:
                    my_gmice = gmice
            a = {'fgrid': fgrid,
                 'dx': gdict['mean_metadata']['dx'],
                 'dy': gdict['mean_metadata']['dy'],
                 'xmin': gdict['mean_metadata']['xmin'],
                 'ymax': gdict['mean_metadata']['ymax'],
                 'contour_levels': contour_levels,
                 'tdir': tdir,
                 'fname': fname,
                 'schema': schema,
                 'imt': imt,
                 'gmice': my_gmice,
                 'gdict': gdict}
            alist.append(a)
            copyfile(os.path.join(smdata, 'WGS1984.prj'),
                     os.path.join(tdir, fname + '.prj'))
            lyrfile = os.path.join(smdata, fname + '.lyr')
            if not os.path.isfile(lyrfile):
                logger.warning("No " + fname + ".lyr file in " + smdata)
            else:
                copyfile(lyrfile, os.path.join(tdir, fname + '.lyr'))
            xmlfile = os.path.join(smdata, fname + '.shp.xml')
            if not os.path.isfile(xmlfile):
                logger.warning("No " + fname + ".shp.xml file in " + smdata)
            else:
                copyfile(xmlfile, os.path.join(tdir, fname + '.shp.xml'))

        worker = partial(make_shape_files, method=method)

        if max_workers > 0:
            with cf.ProcessPoolExecutor(max_workers=max_workers) as ex:
                results = ex.map(worker, alist)
                list(results)
        else:
            for adict in alist:
                worker(adict)

        zfilename = os.path.join(datadir, 'shape.zip')
        zfile = zipfile.ZipFile(zfilename, mode='w',
                                compression=zipfile.ZIP_DEFLATED)
        filelist = []
        for (dirpath, dirnames, filenames) in os.walk(tdir):
            filelist.extend(filenames)
            break
        for sfile in filelist:
            zfile.write(os.path.join(tdir, sfile), sfile)
        zfile.close()


def make_shape_files(adict, method='pcontour'):
    fgrid = adict['fgrid']
    dx = adict['dx']
    dy = adict['dy']
    xmin = adict['xmin']
    ymax = adict['ymax']
    contour_levels = adict['contour_levels']
    tdir = adict['tdir']
    fname = adict['fname']
    schema = adict['schema']
    gdict = adict['gdict']
    imt = adict['imt']
    gmice = adict['gmice']

    if method == 'pcontour':
        gjson = pcontour(fgrid, dx, dy, xmin, ymax, contour_levels, 4, 0, fmt=1)
        features = gjson['features']
    elif method == 'skimage':
        features = contour(gdict, imt, 10, gmice)
    else:
        raise ValueError('Unknown contour method.')
    with fiona.open(os.path.join(tdir, fname + '.shp'), 'w',
                    'ESRI Shapefile', schema) as c:
        for jobj in features:
            c.write(jobj)
