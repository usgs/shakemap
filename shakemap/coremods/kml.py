# stdlib imports
import os
import os.path
import zipfile
import shutil
import re

# third party imports
from PIL import Image
from lxml import etree
import numpy as np
from scipy.ndimage.filters import median_filter
import simplekml as skml
import fiona
import cartopy.io.shapereader as shpreader
from shapely.geometry import shape

# local imports
from mapio.geodict import GeoDict
from impactutils.io.smcontainers import ShakeMapOutputContainer
from .base import CoreModule, Contents
from shakemap.utils.config import (get_config_paths)
from shakelib.plotting.contour import contour
from impactutils.colors.cpalette import ColorPalette
from mapio.grid2d import Grid2D
from shakemap.c.pcontour import pcontour

OVERLAY_IMG = 'ii_overlay.png'
OVERLAY_KML = 'overlay.kml'
STATION_KML = 'stations.kml'
CONTOUR_KML = 'mmi_contour.kml'
POLYGON_KML = 'polygons_mi.kml'
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

color_hash = {
    "0": "ffffffff",
    "0.0": "ffffffff",
    "0.1": "ffffffff",
    "0.2": "ffffffff",
    "0.3": "ffffffff",
    "0.4": "ffffffff",
    "0.5": "ffffffff",
    "0.6": "ffffffff",
    "0.7": "ffffffff",
    "0.8": "ffffffff",
    "0.9": "ffffffff",
    "1": "ffffffff",
    "1.0": "ffffffff",
    "1.1": "fffffaf9",
    "1.2": "fffff5f2",
    "1.3": "fffff0ec",
    "1.4": "ffffebe5",
    "1.5": "ffffe5df",
    "1.6": "ffffe0d9",
    "1.7": "ffffdbd2",
    "1.8": "ffffd6cc",
    "1.9": "ffffd1c5",
    "2": "ffffccbf",
    "2.0": "ffffccbf",
    "2.1": "ffffcfbc",
    "2.2": "ffffd1b9",
    "2.3": "ffffd4b6",
    "2.4": "ffffd6b3",
    "2.5": "ffffd9b0",
    "2.6": "ffffdcac",
    "2.7": "ffffdea9",
    "2.8": "ffffe1a6",
    "2.9": "ffffe3a3",
    "3": "ffffe6a0",
    "3.0": "ffffe6a0",
    "3.1": "ffffe99d",
    "3.2": "ffffeb9a",
    "3.3": "ffffee96",
    "3.4": "fffff093",
    "3.5": "fffff390",
    "3.6": "fffff58d",
    "3.7": "fffff88a",
    "3.8": "fffffa86",
    "3.9": "fffffd83",
    "4": "ffffff80",
    "4.0": "ffffff80",
    "4.1": "fff4ff7f",
    "4.2": "ffe9ff7f",
    "4.3": "ffdfff7e",
    "4.4": "ffd4ff7e",
    "4.5": "ffc9ff7d",
    "4.6": "ffbeff7c",
    "4.7": "ffb3ff7c",
    "4.8": "ffa9ff7b",
    "4.9": "ff9eff7b",
    "5": "ff93ff7a",
    "5.0": "ff93ff7a",
    "5.1": "ff84ff87",
    "5.2": "ff76ff95",
    "5.3": "ff67ffa2",
    "5.4": "ff58ffaf",
    "5.5": "ff4affbc",
    "5.6": "ff3bffca",
    "5.7": "ff2cffd7",
    "5.8": "ff1dffe4",
    "5.9": "ff0ffff2",
    "6": "ff00ffff",
    "6.0": "ff00ffff",
    "6.1": "ff00faff",
    "6.2": "ff00f4ff",
    "6.3": "ff00efff",
    "6.4": "ff00e9ff",
    "6.5": "ff00e4ff",
    "6.6": "ff00deff",
    "6.7": "ff00d9ff",
    "6.8": "ff00d3ff",
    "6.9": "ff00ceff",
    "7": "ff00c8ff",
    "7.0": "ff00c8ff",
    "7.1": "ff00c3ff",
    "7.2": "ff00bdff",
    "7.3": "ff00b8ff",
    "7.4": "ff00b2ff",
    "7.5": "ff00adff",
    "7.6": "ff00a7ff",
    "7.7": "ff00a2ff",
    "7.8": "ff009cff",
    "7.9": "ff0097ff",
    "8": "ff0091ff",
    "8.0": "ff0091ff",
    "8.1": "ff0083ff",
    "8.2": "ff0074ff",
    "8.3": "ff0066ff",
    "8.4": "ff0057ff",
    "8.5": "ff0049ff",
    "8.6": "ff003aff",
    "8.7": "ff002cff",
    "8.8": "ff001dff",
    "8.9": "ff000fff",
    "9": "ff0000ff",
    "9.0": "ff0000ff",
    "9.1": "ff0000fa",
    "9.2": "ff0000f4",
    "9.3": "ff0000ef",
    "9.4": "ff0000e9",
    "9.5": "ff0000e4",
    "9.6": "ff0000de",
    "9.7": "ff0000d9",
    "9.8": "ff0000d3",
    "9.9": "ff0000ce",
    "10": "ff0000c8",
    "10.0": "ff0000c8",
    "10.1": "ff0000c6",
    "10.2": "ff0000c3",
    "10.3": "ff0000c1",
    "10.4": "ff0000be",
    "10.5": "ff0000bc",
    "10.6": "ff0000ba",
    "10.7": "ff0000b7",
    "10.8": "ff0000b5",
    "10.9": "ff0000b2",
    "11": "ff0000b0",
    "11.0": "ff0000b0"}

arabic2roman = {
    "1": "I",
    "2": "II",
    "3": "III",
    "4": "IV",
    "5": "V",
    "6": "VI",
    "7": "VII",
    "8": "VIII",
    "9": "IX",
    "10": "X"}


class KMLModule(CoreModule):
    """
    kml -- Generate KML/KMZ files for ShakeMap.
    """

    command_name = 'kml'
    targets = [r'products/shakemap\.kmz']
    dependencies = [('products/shake_result.hdf', True)]

    def __init__(self, eventid):
        """
        Instantiate a KMLModule class with an event ID.
        """
        super(KMLModule, self).__init__(eventid)
        self.contents = Contents('Ground MOtion KMZ File', 'kml', eventid)

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

        # call create_kmz function
        create_kmz(container, datadir, self.logger, self.contents)

        container.close()


def create_kmz(container, datadir, logger, contents):
    # we're going to combine all these layers into one KMZ file.
    kmz_contents = []

    # create the kml text
    info = container.getMetadata()
    eid = info['input']['event_information']['event_id']
    mag = info['input']['event_information']['magnitude']
    timestr = info['input']['event_information']['origin_time']
    namestr = 'ShakeMap %s M%s %s' % (eid, mag, timestr)
    document = skml.Kml(name=namestr)
    nlc = skml.NetworkLinkControl(minrefreshperiod=300)
    document.networklinkcontrol = nlc

    set_look(document, container)

    # create intensity overlay
    logger.debug('Creating intensity overlay...')
    overlay_image = create_overlay(container, datadir, document)
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

    # create MMI polygon kml
    logger.debug('Creating polygon KML...')
    create_polygons(container, document)
    logger.debug('Created polygon KML')

    # create epicenter KML
    logger.debug('Creating epicenter KML...')
    create_epicenter(container, document)
    logger.debug('Created epicenter KML')

    # place ShakeMap legend on the screen
    legend_file = place_legend(datadir, document)
    kmz_contents.append(legend_file)

    # Write the uber-kml file
    kmlfile = os.path.join(datadir, KML_FILE)
    document.save(kmlfile)
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

    ftype = 'application/vnd.google-earth.kml+xml'
    contents.addFile('shakemap_kmz', 'ShakeMap Overview KMZ',
                     'ShakeMap Overview KMZ.',
                     'shakemap.kmz', ftype)

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
    icon = skml.Icon(href=LEGEND)
    overlayxy = skml.OverlayXY(x=0, y=90, xunits='pixels', yunits='pixels')
    screenxy = skml.ScreenXY(x=5, y=1, xunits='pixels', yunits='fraction')
    size = skml.Size(x=0, y=0, xunits='pixels', yunits='pixels')
    document.newscreenoverlay(name='Intensity Legend', icon=icon,
                              overlayxy=overlayxy,
                              screenxy=screenxy,
                              size=size)

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
    icon = skml.Icon(href=EPICENTER_URL)
    iconstyle = skml.IconStyle(icon=icon)
    style = skml.Style(iconstyle=iconstyle)

    info = container.getMetadata()
    lon = info['input']['event_information']['longitude']
    lat = info['input']['event_information']['latitude']
    point = document.newpoint(name='Earthquake Epicenter',
                              coords=[(lon, lat)],
                              visibility=0)
    point.style = style


def create_polygons(container, document):

    component = container.getComponents('MMI')[0]
    gdict = container.getIMTGrids("MMI", component)
    fgrid = median_filter(gdict['mean'], size=10)
    cont_min = np.floor(np.min(fgrid)) - 0.5
    if cont_min < 0:
        cont_min = 0.5
    cont_max = np.ceil(np.max(fgrid)) + 0.5
    if cont_max > 10.5:
        cont_max = 10.5
    contour_levels = np.arange(cont_min, cont_max, 1, dtype=np.double)
    gjson = pcontour(fgrid,
                     gdict['mean_metadata']['dx'],
                     gdict['mean_metadata']['dy'],
                     gdict['mean_metadata']['xmin'],
                     gdict['mean_metadata']['ymax'],
                     contour_levels, 4, 0)

    folder = document.newfolder(name="MMI Polygons")

    for feature in gjson['features']:
        cv = feature['properties']['value']
        f = folder.newfolder(name="MMI %g Polygons" % cv)
        color = color_hash["%g" % cv]
        name = "MMI %g Polygon" % cv
        s = skml.PolyStyle(fill=1, outline=0, color=color,
                           colormode='normal')
        for plist in feature['geometry']['coordinates']:
            ib = []
            for i, coords in enumerate(plist):
                if i == 0:
                    ob = coords
                else:
                    ib.append(coords)
            p = f.newpolygon(outerboundaryis=ob, innerboundaryis=ib,
                             name=name, visibility=0)
            p.style.polystyle = s

    # Make the polygon labels
    cont_min = np.floor(np.min(fgrid))
    cont_max = np.ceil(np.max(fgrid))
    contour_levels = np.arange(cont_min, cont_max, 1, dtype=np.double)
    gjson = pcontour(fgrid,
                     gdict['mean_metadata']['dx'],
                     gdict['mean_metadata']['dy'],
                     gdict['mean_metadata']['xmin'],
                     gdict['mean_metadata']['ymax'],
                     contour_levels, 2, 0)

    f = folder.newfolder(name="MMI Labels")
    ic = skml.IconStyle(scale=0)
    for feature in gjson['features']:
        cv = "%g" % feature['properties']['value']
        if cv.endswith(".5"):
            continue
        for coords in feature['geometry']['coordinates']:
            lc = len(coords)
            if lc < 150:
                continue
            if coords[0][0] == coords[-1][0] and coords[0][1] == coords[-1][1]:
                if lc < 500:
                    dopts = [0, int(lc/2)]
                elif lc < 1000:
                    dopts = [0, int(lc/3), int(2*lc/3)]
                else:
                    dopts = [0, int(lc/4), int(lc/2), int(3*lc/4)]
            else:
                dopts = [int(lc/2)]
            for i in dopts:
                p = f.newpoint(name=arabic2roman[cv], coords=[coords[i]],
                               visibility=0)
                p.style.iconstyle = ic


def create_contours(container, document):
    """Create a KML file containing MMI contour lines.

    Args:
        container (ShakeMapOutputContainer): Results of model.conf.
        datadir (str): Path to data directory where output KMZ will be written.
        document (Element): LXML KML Document element.
    """
    # TODO - label contours? gx:labelVisibility doesn't seem to be working...

    folder = document.newfolder(name='Contours', visibility=0)
    mmi_line_styles = create_line_styles()
    pgm_line_style = skml.Style(linestyle=skml.LineStyle(width=3))
    ic = skml.IconStyle(scale=0)

    component = list(container.getComponents())[0]
    imts = container.getIMTs(component)
    for imt in imts:
        line_strings = contour(container.getIMTGrids(imt, component), imt,
                               DEFAULT_FILTER_SIZE, None)
        # make a folder for the contours
        imt_folder = folder.newfolder(name='%s Contours' % imt,
                                      visibility=0)

        for line_string in line_strings:
            if imt == 'MMI':
                val = '%.1f' % line_string['properties']['value']
            else:
                val = '%g' % line_string['properties']['value']
            line_list = []
            for segment in line_string['geometry']['coordinates']:
                ctext = []
                for vertex in segment:
                    ctext.append((vertex[0], vertex[1]))
                ls = skml.LineString(coords=ctext)
                line_list.append(ls)
                lc = len(ctext)
                if lc < 10:
                    dopts = []
                elif (ctext[0][0] == ctext[-1][0] and
                      ctext[0][1] == ctext[-1][1]):
                    if lc < 30:
                        dopts = [0, int(lc/2)]
                    elif lc < 60:
                        dopts = [0, int(lc/3), int(2*lc/3)]
                    else:
                        dopts = [0, int(lc/4), int(lc/2), int(3*lc/4)]
                else:
                    dopts = [int(lc/2)]
                for i in dopts:
                    p = imt_folder.newpoint(name=val, coords=[ctext[i]],
                                            visibility=0)
                    p.style.iconstyle = ic
            mg = imt_folder.newmultigeometry(geometries=line_list,
                                             visibility=0,
                                             name="%s %s" % (imt, val))
            if imt == 'MMI':
                mg.style = mmi_line_styles[val]
            else:
                mg.style = pgm_line_style


def set_look(document, container):
    """Set the view location, altitude, and angle.

    Args:
        document (Element): LXML KML Document element.
        container (ShakeMapOutputContainer): Results of model.conf.

    """
    # set the view so that we're looking straight down
    info = container.getMetadata()
    lon = info['input']['event_information']['longitude']
    lat = info['input']['event_information']['latitude']
    document.lookat = skml.LookAt(longitude=lon, latitude=lat,
                                  altitude='%i' % LOOKAT_ALTITUDE,
                                  altitudemode='absolute',
                                  tilt=0, heading=0)


def create_line_styles():
    """Create line styles for contour KML.

    Args:
    """
    line_styles = {}
    cpalette = ColorPalette.fromPreset('mmi')
    for mmi in np.arange(0, 11, 0.5):
        pid = '%.1f' % mmi
        rgb = cpalette.getDataColor(mmi, color_format='hex')
        line_style = skml.LineStyle(color=flip_rgb(rgb),
                                    width=2.0)
        style = skml.Style(linestyle=line_style)
        line_styles[pid] = style
    return line_styles


def create_overlay(container, datadir, document):
    """Create a KML file and intensity map.

    Args:
        container (ShakeMapOutputContainer): Results of model.conf.
        datadir (str): Path to data directory where output KMZ will be written.
        document (SubElement): KML document where the overlay tags should go.
    Returns:
        tuple: (Path to output KMZ file, Path to output overlay image)
    """
    # create the overlay image file
    overlay_img_file = os.path.join(datadir, OVERLAY_IMG)
    geodict = create_overlay_image(container, overlay_img_file)
    box = skml.LatLonBox(north=geodict.ymax, south=geodict.ymin,
                         east=geodict.xmax, west=geodict.xmin)
    icon = skml.Icon(refreshinterval=300,
                     refreshmode='onInterval',
                     href=OVERLAY_IMG)
    document.newgroundoverlay(name='IntensityOverlay',
                              color='ffffffff',
                              draworder=0,
                              latlonbox=box,
                              icon=icon)

    return overlay_img_file


def create_overlay_image(container, filename):
    """Create a semi-transparent PNG image of intensity.

    Args:
        container (ShakeMapOutputContainer): Results of model.conf.
        filename (str): Path to desired output PNG file.
    Returns:
        GeoDict: GeoDict object for the intensity grid.
    """
    # extract the intensity data from the container
    comp = container.getComponents('MMI')[0]
    imtdict = container.getIMTGrids('MMI', comp)
    mmigrid = imtdict['mean']
    gd = GeoDict(imtdict['mean_metadata'])
    imtdata = mmigrid.copy()
    rows, cols = imtdata.shape

    # get the intensity colormap
    palette = ColorPalette.fromPreset('mmi')

    # map intensity values into
    # RGBA array
    rgba = palette.getDataColor(imtdata, color_format='array')

    # set the alpha value to 255 wherever we have MMI 0
    rgba[imtdata <= 1.5] = 0

    if 'CALLED_FROM_PYTEST' not in os.environ:
        # mask off the areas covered by ocean
        oceans = shpreader.natural_earth(category='physical',
                                         name='ocean',
                                         resolution='10m')
        bbox = (gd.xmin, gd.ymin, gd.xmax, gd.ymax)
        with fiona.open(oceans) as c:
            tshapes = list(c.items(bbox=bbox))
            shapes = []
            for tshp in tshapes:
                shapes.append(shape(tshp[1]['geometry']))
            if len(shapes):
                oceangrid = Grid2D.rasterizeFromGeometry(shapes, gd,
                                                         fillValue=0.0)
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

    # get a color palette object to convert intensity values to
    # html colors
    cpalette = ColorPalette.fromPreset('mmi')

    # get the station data from the container
    station_dict = container.getStationDict()

    # Group the MMI and instrumented stations separately
    mmi_folder = document.newfolder(name="Macroseismic Stations",
                                    visibility=0)
    ins_folder = document.newfolder(name="Instrumented Stations",
                                    visibility=0)

    for station in station_dict['features']:
        intensity = get_intensity(station)
        rgb = cpalette.getDataColor(intensity, color_format='hex')
        color = flip_rgb(rgb)
        if station['properties']['station_type'] == 'seismic':
            style_map = create_styles(document, TRIANGLE, 0.6, 0.8, color)
            make_placemark(ins_folder, station, cpalette, style_map)
        else:
            style_map = create_styles(document, CIRCLE, 0.4, 0.6, color)
            make_placemark(mmi_folder, station, cpalette, style_map)

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


def make_placemark(folder, station, cpalette, style_map):
    """Create a placemark element in station KML.

    Args:
        folder (Element): KML Folder element.
        station (dict): Dictionary containing station data.
        cpalette (ColorPalette): Object allowing user to convert MMI to color.
        style_map (skml.StyleMap): The style map for the station type.
    """
    point = folder.newpoint(name=station['id'],
                            visibility=0,
                            description=get_description_table(station),
                            altitudemode='clampToGround',
                            coords=[tuple(station['geometry']['coordinates'])])
    point.stylemap = style_map


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
    period = re.search(r"\d+\.\d+", imt).group()  # noqa
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


def create_styles(document, icon_text, scale_normal, scale_highlight, color):
    """Create styles/style maps for station KML.

    Args:
        document (Element): LXML KML Document element.
    """
    style_normal = add_icon_style(document, icon_text, scale_normal, 0.0,
                                  color)
    style_highlight = add_icon_style(document, icon_text, scale_highlight, 1.0,
                                     color)

    style_map = skml.StyleMap(normalstyle=style_normal,
                              highlightstyle=style_highlight)
    return style_map


def add_icon_style(document, icon_text, icon_scale, label_scale, color):
    """Create Style tag around Icon in KML.

    Args:
        document (Element): LXML KML Document element.
        icon_text (str): The name of the icon file.
        icon_scale (float): The icon scale.
        label_scale (float): The label scale.
    """
    icon = skml.Icon(href=icon_text)
    icon_style = skml.IconStyle(scale="%.1f" % icon_scale,
                                color=color,
                                icon=icon)
    label_style = skml.LabelStyle(scale='%.1f' % label_scale)
#    list_style = skml.ListStyle(listitemtype='checkHideChildren')
    balloon_style = skml.BalloonStyle(text='$[description]')

#    style = skml.Style(iconstyle=icon_style, labelstyle=label_style,
#                       liststyle=list_style, balloonstyle=balloon_style)
    style = skml.Style(iconstyle=icon_style, labelstyle=label_style,
                       balloonstyle=balloon_style)
    return style


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
