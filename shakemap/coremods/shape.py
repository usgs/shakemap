# stdlib imports
import os.path
from collections import OrderedDict
import logging
import zipfile
import tempfile
from shutil import copyfile

# third party imports
import fiona
from impactutils.io.smcontainers import ShakeMapOutputContainer
import numpy as np

# local imports
from .base import CoreModule
from shakemap.utils.config import (get_data_path,
                                   get_config_paths)
from shakemap.utils.logging import get_logging_config
from shakemap.c.pcontour import pcontour
from shakelib.utils.imt_string import oq_to_file


class ShapeModule(CoreModule):
    """
    shape -- Generate shape files of the ground motion parameters
    """

    command_name = 'shape'
    targets = [r'products/shape\.zip']
    dependencies = [('products/shake_result.hdf', True)]

    # supply here a data structure with information about files that
    # can be created by this module.
    shape_page = {'title': 'GIS Shape Files', 'slug': 'shape files'}
    contents = OrderedDict.fromkeys(['shakemap_shapefiles'])
    ftype = 'application/zip'
    contents['shakemap_shapefiles'] = {'title': 'ShakeMap Shape Files',
                                       'caption': 'Shape Files.',
                                       'page': shape_page,
                                       'formats': [{'filename': 'shape.zip',
                                                    'type': ftype}]
                                       }

    def __init__(self, eventid):
        """
        Instantiate a ShapeModule class with an event ID.
        """
        self._eventid = eventid
        log_config = get_logging_config()
        log_name = log_config['loggers'].keys()[0]
        self.logger = logging.getLogger(log_name)

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

        create_polygons(container, datadir, self.logger)

        container.close()


def create_polygons(container, datadir, logger):

    component = list(container.getComponents())[0]
    imts = container.getIMTs(component)

    schema = {'properties': OrderedDict([('AREA', 'float:13.3'),
                                         ('PERIMETER', 'float:14.3'),
                                         ('PGAPOL_', 'int:12'),
                                         ('PGAPOL_ID', 'int:12'),
                                         ('GRID_CODE', 'int:12'),
                                         ('PARAMVALUE', 'float:14.4')]),
              'geometry': 'Polygon'}

    smdata = os.path.join(get_data_path(), 'gis')
    # Make a directory for the files to live in prior to being zipped
    with tempfile.TemporaryDirectory(dir=datadir) as tdir:
        for imt in imts:
            gdict = container.getIMTGrids(imt, component)
            fgrid = gdict['mean']
            if imt == 'MMI':
                contour_levels = np.arange(0.1, 10.2, 0.2)
                fname = 'mi'
            elif imt == 'PGV':
                fgrid = np.exp(fgrid)
                cont_max = np.ceil(np.max(fgrid)) + 2.0
                contour_levels = np.arange(1.0, cont_max, 2.0)
                fname = 'pgv'
            else:
                fgrid = np.exp(fgrid)
                cont_max = (np.ceil(100 * np.max(fgrid)) + 2.0) / 100.0
                contour_levels = np.arange(0.01, cont_max, 0.02)
                fname = oq_to_file(imt)
            gjson = pcontour(fgrid,
                             gdict['mean_metadata']['dx'],
                             gdict['mean_metadata']['dy'],
                             gdict['mean_metadata']['xmin'],
                             gdict['mean_metadata']['ymax'],
                             contour_levels, 4, 0, fmt=1)
            with fiona.open(os.path.join(tdir, fname + '.shp'), 'w',
                            'ESRI Shapefile', schema) as c:
                for jobj in gjson['features']:
                    c.write(jobj)
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
