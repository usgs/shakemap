#stdlib imports
import sys
import os.path
import json
import logging
import glob
import warnings
import io

#third party imports
import numpy as np
from shakelib.utils.containers import OutputContainer
from mapio.gdal import GDALGrid
from configobj import ConfigObj

#local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths,get_logging_config
from shakelib.utils.imt_string import oq_to_file, file_to_oq

FORMATS = {'shapefile':('ESRI Shapefile','shp'),
           'geojson':('GeoJSON','json')}

DEFAULT_FILTER_SIZE = 10

class RasterModule(CoreModule):
    """
    raster - Generate GIS raster files of all configured IMT values.
    """
    command_name = 'raster'
    def execute(self):
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
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

        # create GIS-readable .flt files of imt and uncertainty
        self.logger.info('Creating GIS grids...')
        layers = config['products']['gisgrids']['layers']
        for layer in layers:
            fileimt = oq_to_file(layer)
            imtdict = container.getIMT(layer,'Larger')
            mean_grid = imtdict['mean']
            std_grid = imtdict['std']
            mean_gdal = GDALGrid.copyFromGrid(mean_grid)
            std_gdal = GDALGrid.copyFromGrid(std_grid)
            mean_fname = os.path.join(datadir,'%s_mean.flt' % fileimt)
            std_fname = os.path.join(datadir,'%s_std.flt' % fileimt)
            self.logger.info('Saving %s...' % mean_fname)
            mean_gdal.save(mean_fname)
            self.logger.info('Saving %s...' % std_fname)
            std_gdal.save(std_fname)

