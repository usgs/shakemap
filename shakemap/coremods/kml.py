# stdlib imports
import os.path
from collections import OrderedDict
import logging
import zipfile
import shutil
import re

# third party imports
from shapely.geometry import shape
import fiona
from shakelib.utils.containers import ShakeMapOutputContainer
from PIL import Image
import configobj
from lxml import etree
import numpy as np

# local imports
from .base import CoreModule
from shakemap.utils.config import (get_config_paths,
                                   get_logging_config,
                                   get_configspec,
                                   get_custom_validator,
                                   config_error)
from shakelib.plotting.contour import contour
from impactutils.colors.cpalette import ColorPalette
from mapio.grid2d import Grid2D

OVERLAY_IMG = 'ii_overlay.png'
OVERLAY_KML = 'overlay.kml'
STATION_KML = 'stations.kml'
CONTOUR_KML = 'mmi_contour.kml'
KMZ_FILE = 'shakemap.kmz'
KML_FILE = 'shakemap.kml'
EPICENTER_URL = \
    'http://maps.google.com/mapfiles/kml/shapes/capital_big_highlight.png'
LEGEND = 'intensity_legend.png'

LOOKAT_ALTITUDE = 500000  # meters

TRIANGLE = 'triangle.png'
CIRCLE = 'circle.png'

IMT_UNITS = {'pga': '%g',
             'pgv': 'cm/sec',
             'sa(0.3)': '%g',
             'sa(1.0)': '%g',
             'sa(3.0)': '%g'}

DEFAULT_FILTER_SIZE = 10


class KMLModule(CoreModule):
    """
    kml -- Generate KML/KMZ files for ShakeMap.
    """

    command_name = 'kml'
    targets = [r'products/shakemap\.kmz']
    dependencies = [('products/shake_result.hdf', True)]

    # supply here a data structure with information about files that
    # can be created by this module.
    kml_page = {'title': 'Ground Motion KMZ File', 'slug': 'kml'}
    contents = OrderedDict.fromkeys(['shakemap_kmz'])
    ftype = 'application/vnd.google-earth.kml+xml'
    contents['shakemap_kmz'] = {'title': 'ShakeMap Overview KMZ',
                                'caption': 'ShakeMap Overview.',
                                'page': kml_page,
                                'formats': [{'filename': 'shakemap.kmz',
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
        spec_file = get_configspec('products')
        validator = get_custom_validator()
        pconfig = configobj.ConfigObj(product_config_file,
                                      configspec=spec_file)
        results = pconfig.validate(validator)
        if not isinstance(results, bool) or not results:
            config_error(pconfig, results)
        oceanfile = pconfig['products']['mapping']['layers']['lowres_oceans']

        # call create_kmz function
        create_kmz(container, datadir, oceanfile, self.logger)

        container.close()


def create_kmz(container, datadir, oceanfile, logger):
    # we're going to combine all these layers into one KMZ file.
    kmz_contents = []

    # create the kml text
    root = etree.Element("kml")
    nlink = etree.SubElement(root, "NetworkLinkControl")
    nperiod = etree.SubElement(nlink, "minRefreshPeriod")
    nperiod.text = '300'
    document = etree.SubElement(root, 'Document')
    name = etree.SubElement(document, 'name')
    info = container.getMetadata()
    eid = info['input']['event_information']['event_id']
    mag = info['input']['event_information']['magnitude']
    timestr = info['input']['event_information']['origin_time']
    namestr = 'ShakeMap %s M%s %s' % (eid, mag, timestr)
    name.text = namestr
    set_look(document, container)

    # create intensity overlay
    logger.debug('Creating intensity overlay...')
    overlay_image = create_overlay(container, oceanfile, datadir, document)
    kmz_contents += [overlay_image]
    logger.debug('Created intensity overlay image %s' % overlay_image)

    # create station kml
    logger.debug('Creating station KML...')
    triangle_file, circle_file = create_stations(container,
                                                 datadir,
                                                 document)
    kmz_contents += [triangle_file, circle_file]
    logger.debug('Created station KML')

    # create contour kml
    logger.debug('Creating contour KML...')
    create_contours(container, document)
    logger.debug('Created contour KML')

    # create epicenter KML
    logger.debug('Creating epicenter KML...')
    create_epicenter(container, document)
    logger.debug('Created epicenter KML')

    # place ShakeMap legend on the screen
    legend_file = place_legend(datadir, document)
    kmz_contents.append(legend_file)

    # Write the uber-kml file
    tree = etree.ElementTree(root)
    kmlfile = os.path.join(datadir, KML_FILE)
    tree.write(kmlfile, encoding='utf-8', xml_declaration=True)
    kmz_contents.append(kmlfile)

    # assemble all the pieces into a KMZ file, and delete source files
    # as we go
    kmzfile = os.path.join(datadir, KMZ_FILE)
    kmzip = zipfile.ZipFile(kmzfile, mode='w',
                            compression=zipfile.ZIP_DEFLATED)
    for kfile in kmz_contents:
        _, arcname = os.path.split(kfile)
        kmzip.write(kfile, arcname=arcname)
        os.remove(kfile)
    kmzip.close()

    logger.debug('Wrote KMZ container file %s' % kmzfile)
    return kmzfile


def place_legend(datadir, document):
    """Place the ShakeMap intensity legend in the upper left corner of
    the viewer's map.

    Args:
        datadir (str): Path to data directory where output KMZ will be written.
        document (Element): LXML KML Document element.

    Returns:
        str: Path to output intensity legend file.
    """
    overlay = etree.SubElement(document, 'ScreenOverlay')
    name = etree.SubElement(overlay, 'name')
    name.text = 'Intensity Legend'
    icon = etree.SubElement(overlay, 'Icon')
    href = etree.SubElement(icon, 'href')
    href.text = LEGEND
    _ = etree.SubElement(overlay, 'overlayXY',
                         x="0", y="90",
                         xunits="pixels", yunits="pixels")
    _ = etree.SubElement(overlay, 'screenXY',
                         x="5",
                         y="1",
                         xunits="pixels",
                         yunits="fraction")
    _ = etree.SubElement(overlay, 'size',
                         x="0",
                         y="0",
                         xunits="pixels",
                         yunits="pixels")

    # we need to find the legend file and copy it to
    # the output directory
    this_dir, _ = os.path.split(__file__)
    data_path = os.path.join(this_dir, '..', 'data', 'mapping')
    legend_file = os.path.join(data_path, LEGEND)
    legdest = os.path.join(datadir, LEGEND)
    shutil.copyfile(legend_file, legdest)
    return legdest


def create_epicenter(container, document):
    """Place a star marker at earthquake epicenter.

    Args:
        container (ShakeMapOutputContainer): Results of model.conf.
        document (Element): LXML KML Document element.

    """
    style = etree.SubElement(document, 'Style', id="epicenterIcon")
    iconstyle = etree.SubElement(style, 'IconStyle')
    icon = etree.SubElement(iconstyle, 'Icon')
    href = etree.SubElement(icon, 'href')
    href.text = EPICENTER_URL
    placemark = etree.SubElement(document, 'Placemark')
    name = etree.SubElement(placemark, 'name')
    name.text = 'Earthquake Epicenter'
    point = etree.SubElement(placemark, 'Point')
    coordinates = etree.SubElement(point, 'coordinates')
    info = container.getMetadata()
    lon = info['input']['event_information']['longitude']
    lat = info['input']['event_information']['latitude']
    coordinates.text = '%s,%s' % (lon, lat)
    styleurl = etree.SubElement(placemark, 'styleUrl')
    styleurl.text = '#epicenterIcon'
    visibility = etree.SubElement(placemark, 'visibility')
    visibility.text = '0'


def create_contours(container, document):
    """Create a KML file containing MMI contour lines.

    Args:
        container (ShakeMapOutputContainer): Results of model.conf.
        datadir (str): Path to data directory where output KMZ will be written.
        document (Element): LXML KML Document element.
    """
    # TODO - label contours? gx:labelVisibility doesn't seem to be working...
    imts = container.getIMTs()
    if 'MMI' not in imts:
        return
    component = container.getComponents('MMI')[0]
    line_strings = contour(container, 'MMI', component, DEFAULT_FILTER_SIZE,
                           None)

    # make a folder for the contours
    folder = etree.SubElement(document, 'Folder')
    name = etree.SubElement(folder, 'name')
    name.text = 'Contours'
    visibility = etree.SubElement(folder, 'visibility')
    visibility.text = '0'

    create_line_styles(folder)
    for line_string in line_strings:
        placemark = etree.SubElement(folder, 'Placemark')
        visibility = etree.SubElement(placemark, 'visibility')
        visibility.text = '0'
        styleurl = etree.SubElement(placemark, 'styleUrl')
        mmi = line_string['properties']['value']
        styleurl.text = '#style_mi_%.1f' % mmi
        name = etree.SubElement(placemark, 'name')
        name.text = 'MMI %.1f' % (mmi)
        geometry = etree.SubElement(placemark, 'MultiGeometry')
        for segment in line_string['geometry']['coordinates']:
            linestring = etree.SubElement(geometry, 'LineString')
            coordinates = etree.SubElement(linestring, 'coordinates')
            ctext = ''
            for vertex in segment:
                lon = vertex[0]
                lat = vertex[1]
                ctext += '%.4f,%.4f,0\n' % (lon, lat)
            coordinates.text = ctext


def set_look(document, container):
    """Set the view location, altitude, and angle.

    Args:
        document (Element): LXML KML Document element.
        container (ShakeMapOutputContainer): Results of model.conf.

    """
    # set the view so that we're looking straight down
    lookat = etree.SubElement(document, 'LookAt')
    info = container.getMetadata()
    lon = info['input']['event_information']['longitude']
    lat = info['input']['event_information']['latitude']
    longitude = etree.SubElement(lookat, 'longitude')
    longitude.text = lon
    latitude = etree.SubElement(lookat, 'latitude')
    latitude.text = lat
    altitude = etree.SubElement(lookat, 'altitude')
    altitude.text = '%i' % LOOKAT_ALTITUDE
    tilt = etree.SubElement(lookat, 'tilt')
    tilt.text = '0'
    altmode = etree.SubElement(lookat, 'altitudeMode')
    altmode.text = 'absolute'


def create_line_styles(document):
    """Create line styles for contour KML.

    Args:
        document (Element): LXML KML Document element.
    """
    gxns = 'http://www.google.com/kml/ext/2.2'
    nsmap = {'gx': gxns}
    cpalette = ColorPalette.fromPreset('mmi')
    for mmi in np.arange(0, 11, 0.5):
        pid = 'style_mi_%.1f' % mmi
        style = etree.SubElement(document, 'Style', id=pid)
        linestyle = etree.SubElement(style, 'LineStyle')
        color = etree.SubElement(linestyle, 'color')
        rgb = cpalette.getDataColor(mmi, color_format='hex')
        color.text = flip_rgb(rgb)
        width = etree.SubElement(linestyle, 'width')
        width.text = '2.0'
        # TODO: this doesn't work!
        vis = etree.SubElement(linestyle,
                               '{%s}labelVisibility' % gxns,
                               nsmap=nsmap)
        vis.text = '1'


def create_overlay(container, oceanfile, datadir, document):
    """Create a KML file and intensity map.

    Args:
        container (ShakeMapOutputContainer): Results of model.conf.
        oceanfile (str): Path to shapefile containing ocean polygons.
        datadir (str): Path to data directory where output KMZ will be written.
        document (SubElement): KML document where the overlay tags should go.
    Returns:
        tuple: (Path to output KMZ file, Path to output overlay image)
    """
    # create the overlay image file
    overlay_img_file = os.path.join(datadir, OVERLAY_IMG)
    geodict = create_overlay_image(container, oceanfile, overlay_img_file)

    overlay = etree.SubElement(document, 'GroundOverlay')
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
    href.text = OVERLAY_IMG
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

    return overlay_img_file


def create_overlay_image(container, oceanfile, filename):
    """Create a semi-transparent PNG image of intensity.

    Args:
        container (ShakeMapOutputContainer): Results of model.conf.
        oceanfile (str): Path to shapefile containing ocean polygons.
        filename (str): Path to desired output PNG file.
    Returns:
        GeoDict: GeoDict object for the intensity grid.
    """
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
        if len(shapes):
            oceangrid = Grid2D.rasterizeFromGeometry(shapes, gd, fillValue=0.0)
            rgba[oceangrid.getData() == 1] = 0

    # save rgba image as png
    img = Image.fromarray(rgba)
    img.save(filename)
    return gd


def create_stations(container, datadir, document):
    """Create a KMZ file containing station KML and necessary icons files.

    Args:
        container (ShakeMapOutputContainer): Results of model.conf.
        datadir (str): Path to data directory where output KMZ will be written.
        document (Element): LXML KML Document element.
    Returns:
        str: Path to output KMZ file.
    """
    create_styles(document)

    # get a color palette object to convert intensity values to
    # html colors
    cpalette = ColorPalette.fromPreset('mmi')

    # get the station data from the container
    station_dict = container.getStationDict()

    # Group the MMI and instrumented stations separately
    mmi_folder = etree.SubElement(document, 'Folder')
    mmi_name = etree.SubElement(mmi_folder, 'name')
    mmi_name.text = 'Macroseismic Stations'
    mmi_vis = etree.SubElement(mmi_folder, 'visibility')
    mmi_vis.text = '0'

    ins_folder = etree.SubElement(document, 'Folder')
    ins_name = etree.SubElement(ins_folder, 'name')
    ins_name.text = 'Instrumented Stations'
    ins_vis = etree.SubElement(ins_folder, 'visibility')
    ins_vis.text = '0'
    for station in station_dict['features']:
        if station['properties']['station_type'] == 'seismic':
            make_placemark(ins_folder, station, cpalette)
        else:
            make_placemark(mmi_folder, station, cpalette)

    # we need to find the triangle and circle icons and copy them to
    # the output directory
    this_dir, _ = os.path.split(__file__)
    data_path = os.path.join(this_dir, '..', 'data', 'mapping')
    triangle_file = os.path.join(data_path, TRIANGLE)
    circle_file = os.path.join(data_path, CIRCLE)
    tridest = os.path.join(datadir, TRIANGLE)
    cirdest = os.path.join(datadir, CIRCLE)
    shutil.copyfile(triangle_file, tridest)
    shutil.copyfile(circle_file, cirdest)

    return (tridest, cirdest)


def make_placemark(document, station, cpalette):
    """Create a placemark element in station KML.

    Args:
        document (Element): LXML KML Document element.
        station (dict): Dictionary containing station data.
        cpalette (ColorPalette): Object allowing user to convert MMI to color.
    """
    placemark = etree.SubElement(document, 'Placemark')
    style_url = etree.SubElement(placemark, 'styleUrl')
    if station['properties']['network'] in ['INTENSITY', 'DYFI', 'CIIM']:
        style_url.text = '#dyfiIconMap'
    else:
        style_url.text = '#stationIconMap'
    style = etree.SubElement(placemark, 'Style')
    icon_style = etree.SubElement(style, 'IconStyle')
    color = etree.SubElement(icon_style, 'color')
    intensity = get_intensity(station)
    rgb = cpalette.getDataColor(intensity, color_format='hex')
    color.text = flip_rgb(rgb)
    name = etree.SubElement(placemark, 'name')
    name.text = station['id']
    visibility = etree.SubElement(placemark, 'visibility')
    visibility.text = '0'
    description = etree.SubElement(placemark, 'description')
    description.text = etree.CDATA(get_description_table(station))
    point = etree.SubElement(placemark, 'Point')
    alt_model = etree.SubElement(point, 'altitudeModel')
    alt_model.text = 'clampToGround'
    coordinates = etree.SubElement(point, 'coordinates')
    coordinates.text = '%.4f,%.4f,0' % (
        tuple(station['geometry']['coordinates']))


def get_intensity(station):
    """Retrieve the intensity value from a station dictionary.

    Args:
        station (dict): Dictionary containing station data.

    Returns:
        float: Intensity value.
    """
    intensity = station['properties']['intensity']
    if intensity == 'null':
        channels = [channel['name']
                    for channel in station['properties']['channels']]
        if 'mmi' not in channels:
            intensity = 0
        else:
            mmi_idx = channels.index('mmi')
            mmid = station['properties']['channels'][mmi_idx]
            intensity = mmid['amplitudes'][0]['value']
    return intensity


def get_description_table(station):
    """Get station description as HTML table.

    Args:
        station (dict): Dictionary containing station data.

    Returns:
        str: String containing HTML table describing station.
    """
    table = etree.Element('table', border='1')

    if station['properties']['network'] in ['INTENSITY', 'DYFI', 'CIIM']:
        network = 'DYFI'
    else:
        network = station['properties']['network']
    make_row(table, 'Network', network)
    make_row(table, 'Station ID', station['id'])
    make_row(table, 'Location', station['properties']['location'])
    make_row(table,
             'Lat',
             '%.4f' % station['geometry']['coordinates'][1])
    make_row(table,
             'Lon',
             '%.4f' % station['geometry']['coordinates'][0])
    make_row(table,
             'Distance to source',
             '%.2f' % station['properties']['distance'])

    intensity_string = '%.1f' % get_intensity(station)
    make_row(table, 'Intensity', intensity_string)

    # get the list of IMTs in the station tag
    imt_list = []
    for channel in station['properties']['channels']:
        for amplitude in channel['amplitudes']:
            if amplitude['name'] == 'mmi':
                continue
            imt_list.append(amplitude['name'])

    for imt in imt_list:
        imt_str = imt_to_string(imt)
        make_row(table, imt_str, get_imt_text(station, imt))

    desc_xml = etree.tostring(table).decode('utf8')
    return desc_xml


def imt_to_string(imt):
    non_spectrals = {'pga': 'PGA',
                     'pgv': 'PGV'}
    if imt in non_spectrals:
        return non_spectrals[imt]
    period = re.search("\d+\.\d+", imt).group()  # noqa
    imt_string = 'PSA %s sec' % period
    return imt_string


def make_row(table, key, value):
    """Create a row in the description table.

    Args:
        table (Element): LXML Element for table tag.
        key (str): Text for left hand column of row.
        value (str): Text for right hand column of row.

    """
    row = etree.SubElement(table, 'tr')
    row_col1 = etree.SubElement(row, 'td')
    row_col1.text = key
    row_col2 = etree.SubElement(row, 'td')
    row_col2.text = value


def get_description(station):
    """Get station description as HTML definition list.

    Args:
        station (dict): Dictionary containing station data.

    Returns:
        str: String containing HTML definition list describing station.
    """
    dl = etree.Element("dl")
    network_dt = etree.SubElement(dl, 'dt')
    network_dt.text = 'NETWORK:'
    network_dd = etree.SubElement(dl, 'dd')
    if station['properties']['network'] in ['INTENSITY', 'DYFI', 'CIIM']:
        network_dd.text = 'DYFI'
    else:
        network_dd.text = station['properties']['network']

    station_dt = etree.SubElement(dl, 'dt')
    station_dt.text = 'Station ID:'
    station_dd = etree.SubElement(dl, 'dd')
    station_dd.text = station['id']

    location_dt = etree.SubElement(dl, 'dt')
    location_dt.text = 'Location:'
    location_dd = etree.SubElement(dl, 'dd')
    location_dd.text = station['properties']['location']

    lat_dt = etree.SubElement(dl, 'dt')
    lat_dt.text = 'Lat:'
    lat_dd = etree.SubElement(dl, 'dd')
    lat_dd.text = '%.4f' % station['geometry']['coordinates'][1]

    lon_dt = etree.SubElement(dl, 'dt')
    lon_dt.text = 'Lon:'
    lon_dd = etree.SubElement(dl, 'dd')
    lon_dd.text = '%.4f' % station['geometry']['coordinates'][0]

    distance_dt = etree.SubElement(dl, 'dt')
    distance_dt.text = 'Distance to source:'
    distance_dd = etree.SubElement(dl, 'dd')
    distance_dd.text = '%.2f' % station['properties']['distance']

    intensity_dt = etree.SubElement(dl, 'dt')
    intensity_dt.text = 'Intensity:'
    intensity_dd = etree.SubElement(dl, 'dd')
    intensity_dd.text = '%.1f' % get_intensity(station)

    # get the names of all the IMTs
    # TODO: Get list of IMTs dynamically from the data.
    pga_dt = etree.SubElement(dl, 'dt')
    pga_dt.text = 'PGA:'
    pga_dd = etree.SubElement(dl, 'dd')
    pga_dd.text = get_imt_text(station, 'pga')

    pgv_dt = etree.SubElement(dl, 'dt')
    pgv_dt.text = 'PGV:'
    pgv_dd = etree.SubElement(dl, 'dd')
    pgv_dd.text = get_imt_text(station, 'pgv')

    psa03_dt = etree.SubElement(dl, 'dt')
    psa03_dt.text = 'PSA 0.3 sec:'
    psa03_dd = etree.SubElement(dl, 'dd')
    psa03_dd.text = get_imt_text(station, 'sa(0.3)')

    psa10_dt = etree.SubElement(dl, 'dt')
    psa10_dt.text = 'PSA 1.0 sec:'
    psa10_dd = etree.SubElement(dl, 'dd')
    psa10_dd.text = get_imt_text(station, 'sa(1.0)')

    psa30_dt = etree.SubElement(dl, 'dt')
    psa30_dt.text = 'PSA 3.0 sec:'
    psa30_dd = etree.SubElement(dl, 'dd')
    psa30_dd.text = get_imt_text(station, 'sa(3.0)')

    desc_xml = etree.tostring(dl).decode('utf8')
    return desc_xml


def get_imt_text(station, imt):
    """Get a text string describing the value of input IMT.

    Args:
        station (dict): Dictionary containing station data.
        imt (str): IMT string (pga,pgv, sa(0.3), etc.)

    Returns:
        str: IMT text string (i.e., "7.1 cm/sec")
    """
    imt_max = -1
    units = IMT_UNITS[imt]
    for channel in station['properties']['channels']:
        if channel['name'].endswith('Z'):
            continue
        for amplitude in channel['amplitudes']:
            if amplitude['name'] != imt:
                continue
            imt_value = amplitude['value']
            if imt_value == 'null':
                continue

            if imt_value > imt_max:
                imt_max = imt_value

    if imt_max > -1:
        imt_text = '%.1f %s' % (imt_max, units)
    else:
        imt_text = 'nan %s' % (units)
    return imt_text


def create_styles(document):
    """Create styles/style maps for station KML.

    Args:
        document (Element): LXML KML Document element.
    """
    add_icon_style(document, 'stationIconNormal', TRIANGLE, 0.6, 0.0)
    add_icon_style(document, 'stationIconHighlight', TRIANGLE, 0.8, 1.0)
    add_icon_style(document, 'dyfiIconNormal', CIRCLE, 0.4, 0.0)
    add_icon_style(document, 'dyfiIconHighlight', CIRCLE, 0.6, 1.0)

    add_style_map(document, 'station')
    add_style_map(document, 'dyfi')


def add_icon_style(document, icon_id, icon_text, icon_scale, label_scale):
    """Create Style tag around Icon in KML.

    Args:
        document (Element): LXML KML Document element.
        icon_id (str): Id field in Style tag.
        icon_text (str): The name of the icon file.
        icon_scale (float): The icon scale.
        label_scale (float): The label scale.
    """
    normal_station_style = etree.SubElement(document, 'Style', id=icon_id)
    icon_style = etree.SubElement(normal_station_style, 'IconStyle')
    scale1 = etree.SubElement(icon_style, 'scale')
    scale1.text = '%.1f' % icon_scale
    icon = etree.SubElement(icon_style, 'Icon')
    href = etree.SubElement(icon, 'href')
    href.text = icon_text
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
    """Create the StyleMap tags used by the placemark elements.

    Args:
        document (Element): LXML KML Document element.
        style_type (str): "station" or "dyfi".
    """
    style_map = etree.SubElement(document,
                                 'StyleMap',
                                 id='%sIconMap' % style_type)

    pair1 = etree.SubElement(style_map, 'Pair')
    key1 = etree.SubElement(pair1, 'key')
    key1.text = 'normal'
    style_url = etree.SubElement(pair1, 'styleUrl')
    style_url.text = '#%sIconNormal' % style_type

    pair2 = etree.SubElement(style_map, 'Pair')
    key2 = etree.SubElement(pair2, 'key')
    key2.text = 'highlight'
    style_url = etree.SubElement(pair2, 'styleUrl')
    style_url.text = '#%sIconHighlight' % style_type


def flip_rgb(rgb):
    """Reverse order of RGB hex string, prepend 'ff' for transparency.

    Args:
        rgb: RGB hex string (#E1C2D3)
    Returns:
        str: ABGR hex string (#ffd3c2e1).
    """
    # because Google decided that colors in KML should be
    # specified as ABGR instead of RGBA, we have to reverse the
    # sense of our color.
    # so given #E1C2D3, (red E1, green C2, blue D3) convert to #ffd3c2e1
    # where the first hex pair is transparency and the others are in
    # reverse order.
    abgr = rgb.replace('#', '')
    abgr = '#FF' + abgr[4:6] + abgr[2:4] + abgr[0:2]
    abgr = abgr.lower()
    return abgr
