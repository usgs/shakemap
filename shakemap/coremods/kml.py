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
from PIL import Image
import configobj
from lxml import etree

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths, get_logging_config
from shakemap.utils.utils import path_macro_sub
from impactutils.colors.cpalette import ColorPalette
from mapio.gmt import GMTGrid

FORMATS = {
    'geojson': ('GeoJSON', 'json')
}

DEFAULT_FILTER_SIZE = 10


class KMLModule(CoreModule):
    """
    kml -- Generate KML/KMZ files for ShakeMap.
    """

    command_name = 'kml'

    # supply here a data structure with information about files that
    # can be created by this module.
    kml_page = {'title': 'Ground Motion KML Files', 'slug': 'kml'}
    contents = OrderedDict.fromkeys(['intensity_overlay'])
    contents['intensity_overlay'] = {'title': 'Intensity Overlay',
                                     'caption': 'Intensity Overlay.',
                                     'page': kml_page,
                                     'formats': [{'filename': 'overlay.kmz',
                                                  'type': 'application/vnd.google-earth.kml+xml'}
                                                 ]
                                     }

    def __init__(self, eventid, filter=None):
        """
        Instantiate a KMLModule class with an event ID.
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
        Create KML files.

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
            raise NotImplementedError('kml module can only contour '
                                      'gridded data, not sets of points')

        # find the topography grid configured for the system
        product_config_file = os.path.join(
            install_path, 'config', 'products.conf')
        pconfig = configobj.ConfigObj(product_config_file)
        topofile = pconfig['products']['mapping']['layers']['topography']
        topofile = path_macro_sub(topofile, install_path, data_path)

        # create intensity overlay
        self.logger.debug('Creating intensity overlay...')
        overlay_file = create_overlay_kml(container, topofile, datadir)
        self.logger.debug('Wrote intensity overlay file %s' % overlay_file)


def create_overlay_kml(container, topofile, datadir):
    overlay_img_file = os.path.join(datadir, 'ii_overlay.png')
    create_overlay_image(container, topofile, filename)
    root = etree.Element("kml")
    nlink = etree.SubElement(root, "NetworkLinkControl")
    nperiod = etree.SubElement(nlink, "minRefreshPeriod")


def create_overlay_image(container, topofile, filename):
    # extract the intensity data from the container
    comp = container.getComponents('MMI')[0]
    imtdict = container.getIMTGrids('MMI', comp)
    mmigrid = imtdict['mean']
    gd = mmigrid.getGeoDict()
    imtdata = mmigrid.getData().copy()
    rows, cols = imtdata.shape

    # get the intensity colormap
    palette = ColorPalette.fromPreset('mmi')

    # map intensity values into
    # RGBA array
    rgba = palette.getDataColor(imtdata, color_format='array')

    # set the alpha value to 255 wherever we have MMI 0
    rgba[imtdata <= 1.5] = 0

    # we need to mask off the areas covered by water
    # get topo layer and project it
    topogrid = GMTGrid.load(topofile, samplegeodict=gd, resample=True)
    topodata = topogrid.getData().copy()
    rgba[topodata <= 0] = 0

    # save rgba image as png
    img = Image.fromarray(rgba)
    img.save(filename)
