# stdlib imports
import os.path
from collections import OrderedDict
import logging
import time
import zipfile

# third party imports
from shapely.geometry import shape
import fiona
from shakelib.utils.containers import ShakeMapOutputContainer
from PIL import Image
import configobj
from lxml import etree

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths, get_logging_config
from shakemap.utils.utils import path_macro_sub
from impactutils.colors.cpalette import ColorPalette
from mapio.grid2d import Grid2D

IMG_FILE = 'ii_overlay.png'
KML_FILE = 'overlay.kml'


class KMLModule(CoreModule):
    """
    kml -- Generate KML/KMZ files for ShakeMap.
    """

    command_name = 'kml'

    # supply here a data structure with information about files that
    # can be created by this module.
    kml_page = {'title': 'Ground Motion KML Files', 'slug': 'kml'}
    contents = OrderedDict.fromkeys(['intensity_overlay'])
    ftype = 'application/vnd.google-earth.kml+xml'
    contents['intensity_overlay'] = {'title': 'Intensity Overlay',
                                     'caption': 'Intensity Overlay.',
                                     'page': kml_page,
                                     'formats': [{'filename': 'overlay.kmz',
                                                  'type': ftype}
                                                 ]
                                     }

    def __init__(self, eventid):
        """
        Instantiate a KMLModule class with an event ID.
        """
        self._eventid = eventid
        log_config = get_logging_config()
        log_name = log_config['loggers'].keys()[0]
        self.logger = logging.getLogger(log_name)

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

        # find the low res ocean vector dataset
        product_config_file = os.path.join(
            install_path, 'config', 'products.conf')
        pconfig = configobj.ConfigObj(product_config_file)
        oceanfile = pconfig['products']['mapping']['layers']['lowres_oceans']
        oceanfile = path_macro_sub(oceanfile, install_path, data_path)

        # create intensity overlay
        self.logger.debug('Creating intensity overlay...')
        kmzfile = overlay_file = create_overlay_kml(
            container, oceanfile, datadir)
        self.logger.debug('Wrote intensity overlay file %s' % kmzfile)


def create_overlay_kml(container, oceanfile, datadir):
    # create the overlay image file
    overlay_img_file = os.path.join(datadir, IMG_FILE)
    geodict = create_overlay_image(container, oceanfile, overlay_img_file)

    # create the kml text
    root = etree.Element("kml")
    nlink = etree.SubElement(root, "NetworkLinkControl")
    nperiod = etree.SubElement(nlink, "minRefreshPeriod")
    nperiod.text = '300'
    overlay = etree.SubElement(root, 'GroundOverlay')
    name = etree.SubElement(overlay, 'name')
    name.text = 'Intensity Overlay'
    color = etree.SubElement(overlay, 'color')
    color.text = 'ffffffff'
    draw_order = etree.SubElement(overlay, 'drawOrder')
    draw_order.text = '0'
    icon = etree.SubElement(overlay, 'Icon')
    interval = etree.SubElement(icon, 'refreshInterval')
    interval.text = '300'
    mode = etree.SubElement(icon, 'refreshMode')
    mode.text = 'onInterval'
    href = etree.SubElement(icon, 'href')
    href.text = 'ii_overlay.png'
    box = etree.SubElement(overlay, 'LatLonBox')
    xmin, xmax, ymin, ymax = (geodict.xmin, geodict.xmax,
                              geodict.ymin, geodict.ymax)
    north = etree.SubElement(box, 'north')
    north.text = '%.4f' % ymax
    south = etree.SubElement(box, 'south')
    south.text = '%.4f' % ymin
    east = etree.SubElement(box, 'east')
    east.text = '%.4f' % xmax
    west = etree.SubElement(box, 'west')
    west.text = '%.4f' % xmin

    # write the kml file
    kmlfile = os.path.join(datadir, KML_FILE)
    tree = etree.ElementTree(root)
    tree.write(kmlfile, pretty_print=True)

    # create the overlay.kmz file
    kmzfile = os.path.join(datadir, 'overlay.kmz')
    kmz = zipfile.ZipFile(kmzfile, mode='w', compression=zipfile.ZIP_DEFLATED)
    kmz.write(kmlfile, arcname=KML_FILE)
    kmz.write(overlay_img_file, arcname=IMG_FILE)
    kmz.close()
    os.remove(kmlfile)
    os.remove(overlay_img_file)

    return kmzfile


def create_overlay_image(container, oceanfile, filename):
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

    # mask off the areas covered by ocean
    bbox = (gd.xmin, gd.ymin, gd.xmax, gd.ymax)
    with fiona.open(oceanfile) as c:
        tshapes = list(c.items(bbox=bbox))
        shapes = []
        for tshp in tshapes:
            shapes.append(shape(tshp[1]['geometry']))
        oceangrid = Grid2D.rasterizeFromGeometry(shapes, gd, fillValue=0.0)
        rgba[oceangrid.getData() == 1] = 0

    # save rgba image as png
    img = Image.fromarray(rgba)
    img.save(filename)
    return gd


def create_station_kml()


def create_styles(document):
    triangle = 'triangle.png'
    circle = 'icon18.png'
    add_icon_style(document, 'stationIconNormal', triangle, 0.6, 0.0)
    add_icon_style(document, 'stationIconHighLight', triangle, 0.8, 1.0)
    add_icon_style(document, 'dyfiIconNormal', circle, 0.4, 0.0)
    add_icon_style(document, 'stationIconHighLight', circle, 0.6, 1.0)

    add_style_map(document, 'station')
    add_style_map(document, 'dyfi')


def add_icon_style(document, id, icon, icon_scale, label_scale):
    normal_station_style = etree.SubElement(document, 'Style', id=id)
    icon_style = etree.SubElement(normal_station_style, 'IconStyle')
    scale1 = etree.SubElement(icon_style1, 'scale')
    scale1.text = '%.1f' % icon_scale
    icon = etree.SubElement(icon_style, 'Icon')
    icon.text = icon
    label_style = etree.SubElement(normal_station_style, 'LabelStyle')
    scale2 = etree.SubElement(label_style, 'scale')
    scale2.text = '%.1f' % label_scale
    list_style = etree.SubElement(normal_station_style, 'ListStyle')
    item_type = etree.SubElement(list_style, 'listItemType')
    item_type.text = 'checkHideChildren'
    balloon_style = etree.SubElement(normal_station_style, 'BalloonStyle')
    text = etree.SubElement(balloon_style, 'text')
    text.text = '$[description]'


def add_style_map(document, style_type):
    style_map = etree.SubElement(document,
                                 'StyleMap',
                                 id='%sIconMap' % style_type)

    pair1 = etree.SubElement(style_map, 'Pair')
    key1 = etree.SubElement(pair1, 'key')
    key1.text = 'normal'
    style_url = etree.SubElement(pair, 'styleUrl')
    style_url.text = '#%sIconNormal' % style_type

    pair2 = etree.SubElement(style_map, 'Pair')
    key2 = etree.SubElement(pair2, 'key')
    key2.text = 'highlight'
    style_url = etree.SubElement(pair, 'styleUrl')
    style_url.text = '#%sIconHighlight' % style_type


def create_styles(document):
    triangle = 'triangle.png'
    circle = 'icon18.png'
    add_icon_style(document, 'stationIconNormal', triangle, 0.6, 0.0)
    add_icon_style(document, 'stationIconHighLight', triangle, 0.8, 1.0)
    add_icon_style(document, 'dyfiIconNormal', circle, 0.4, 0.0)
    add_icon_style(document, 'stationIconHighLight', circle, 0.6, 1.0)

    add_style_map(document, 'station')
    add_style_map(document, 'dyfi')
