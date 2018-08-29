# stdlib imports
from datetime import datetime
import os.path
import time

# third party imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LightSource

import cartopy.crs as ccrs  # projections
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import ShapelyFeature
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader

from shapely.geometry import shape as sShape
from shapely.geometry import Polygon as sPolygon
from shapely.geometry import GeometryCollection
from shapely.geometry import mapping

import pyproj
import fiona
from openquake.hazardlib import imt

# neic imports
from mapio.gmt import GMTGrid
from impactutils.mapping.mercatormap import MercatorMap
from impactutils.mapping.city import Cities
from impactutils.colors.cpalette import ColorPalette
from impactutils.mapping.scalebar import draw_scale

# local imports
from shakelib.rupture.point_rupture import PointRupture
from shakelib.rupture import constants
from shakelib.plotting.contour import contour
from shakelib.utils.imt_string import oq_to_file
from shakemap.utils.utils import get_object_from_config

# define some constants
WATERCOLOR = '#7AA1DA'
FIGWIDTH = 7.0
FILTER_SMOOTH = 5.0
XOFFSET = 4  # how many pixels between the city dot and the city text
VERT_EXAG = 0.1  # what is the vertical exaggeration for hillshade

# define the zorder values for various map components
# all of the zorder values for different plotted parameters
# elements with a higher zorder will plot on top of lower elements.
IMG_ZORDER = 1
ROAD_ZORDER = 5
COAST_ZORDER = 11
CONTOUR_ZORDER = 800
DASHED_CONTOUR_ZORDER = 1002
OCEAN_ZORDER = 1000
BORDER_ZORDER = 1001
FAULT_ZORDER = 1100
EPICENTER_ZORDER = 1100
STATIONS_ZORDER = 1150
CITIES_ZORDER = 1200
GRATICULE_ZORDER = 1200
SCALE_ZORDER = 1500
SCENARIO_ZORDER = 2000
AXES_ZORDER = 3000

# default font for cities
# DEFAULT_FONT = 'DejaVu Sans'
DEFAULT_FONT = 'Bitstream Vera Sans'

# what fraction of the width/height of the map can the contour labels be
# from the edge before we decide we don't want to draw it there?
MAP_FRAC = 10

# define dictionary of MMI integer values to Roman numeral equivalents
MMI_LABELS = {'1': 'I',
              '2': 'II',
              '3': 'III',
              '4': 'IV',
              '5': 'V',
              '6': 'VI',
              '7': 'VII',
              '8': 'VIII',
              '9': 'IX',
              '10': 'X'}

DEG2KM = 111.191


def get_projected_grids(imtgrid, topobase, projstr):
    # resample shakemap to topogrid
    # get the geodict for the topo file
    topodict = topobase.getGeoDict()

    # get a geodict that is aligned with topo, but inside shakemap
    sampledict = topodict.getBoundsWithin(imtgrid.getGeoDict())

    simtgrid = imtgrid.interpolateToGrid(sampledict)

    # get topo layer and project it
    t1 = time.time()
    topogrid = topobase.interpolateToGrid(sampledict)
    t2 = time.time()
    print('Topo grid took %.2f seconds to load.' % (t2-t1))

    # resampling 32bit floats gives odd results... upcasting to 64bit
    topogrid._data = topogrid._data.astype(np.float64)

    # project the topography data
    ptopogrid = topogrid.project(projstr, method='bilinear')

    # project the MMI data
    pimtgrid = simtgrid.project(projstr)

    return (pimtgrid, ptopogrid)


def get_map_info(gd, center_lat):
    # define the map
    # first cope with stupid 180 meridian
    height = (gd.ymax-gd.ymin)*DEG2KM
    if gd.xmin < gd.xmax:
        width = (gd.xmax-gd.xmin)*np.cos(np.radians(center_lat))*DEG2KM
        xmin, xmax, ymin, ymax = (gd.xmin, gd.xmax, gd.ymin, gd.ymax)
    else:
        xmin, xmax, ymin, ymax = (gd.xmin, gd.xmax, gd.ymin, gd.ymax)
        xmax += 360
        width = ((gd.xmax+360) - gd.xmin) * \
            np.cos(np.radians(center_lat))*DEG2KM

    aspect = width/height

    # if the aspect is not 1, then trim bounds in x or y direction
    # as appropriate
    if width > height:
        dw = (width - height)/2.0  # this is width in km
        xmin = xmin + dw/(np.cos(np.radians(center_lat))*DEG2KM)
        xmax = xmax - dw/(np.cos(np.radians(center_lat))*DEG2KM)
        width = (xmax-xmin)*np.cos(np.radians(center_lat))*DEG2KM
    if height > width:
        dh = (height - width)/2.0  # this is width in km
        ymin = ymin + dh/DEG2KM
        ymax = ymax - dh/DEG2KM
        height = (ymax-ymin)*DEG2KM

    aspect = width/height
    figheight = FIGWIDTH/aspect
    bounds = (xmin, xmax, ymin, ymax)
    figsize = (FIGWIDTH, figheight)
    return (bounds, figsize)


def draw_colorbar(fig, mmimap):
    # making up our own colorbar object here because the default
    # pyplot functionality doesn't seem to do the job.
    cax = fig.add_axes([0.1, 0.035, 0.8, 0.035])
    cax.get_yaxis().set_visible(False)
    cax_xmin, cax_xmax = cax.get_xlim()
    bottom, top = cax.get_ylim()
    plt.xlim(cax_xmin, cax_xmax)
    plt.ylim(bottom, top)
    nsteps = 200
    left = 0
    rights = np.arange(1/nsteps, 1+(1/nsteps), 1/nsteps)
    mmis = np.linspace(1, 10, nsteps)
    for mmi, right in zip(mmis, rights):
        px = [left, right, right, left, left]
        py = [top, top, bottom, bottom, top]
        mmicolor = mmimap.getDataColor(mmi, color_format='hex')
        left = right
        plt.fill(px, py, mmicolor, ec=mmicolor)

    start_loc = (1/nsteps) - (1/nsteps)/2
    end_loc = 1 - (1/nsteps)/2
    locs = np.linspace(start_loc, end_loc, 10)
    labels = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X']

    plt.xticks(locs, labels)


def label_close_to_edge(x, y, xmin, xmax, ymin, ymax):
    # get the map width/heights
    width = xmax - xmin
    height = ymax - ymin

    # start with the first index
    dleft = x - xmin
    dright = xmax - x
    dtop = ymax - y
    dbottom = y - ymin
    return (dleft < width/MAP_FRAC or dright < width/MAP_FRAC or
            dtop < height/MAP_FRAC or dbottom < height/MAP_FRAC)


def draw_graticules(ax, xmin, xmax, ymin, ymax):
    gl = ax.gridlines(draw_labels=True,
                      linewidth=0.5, color='k',
                      alpha=0.5, linestyle='-',
                      zorder=GRATICULE_ZORDER)
    gl.xlabels_top = False
    gl.xlabels_bottom = True
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xlines = True

    # create a dictionary with the intervals we want for a span
    # of degrees.
    spans = {1: 4,
             2: 4,
             3: 2,
             5: 2,
             7: 0.5}

    span_keys = np.array(sorted(list(spans.keys())))
    nearest_span_idx = np.argmin(int((xmax-xmin)) - span_keys)
    interval = spans[span_keys[nearest_span_idx]]

    # let's floor/ceil the edges to nearest 1/interval
    gxmin = np.floor(xmin * interval) / interval
    gxmax = np.ceil(xmax * interval) / interval
    gymin = np.floor(ymin * interval) / interval
    gymax = np.ceil(ymax * interval) / interval

    # check for meridian crossing
    crosses = False
    if gxmax < 0 and gxmax < gxmin:
        crosses = True
        gxmax += 360

    # shakemap way
    ylocs = np.arange(np.floor(gymin), np.ceil(gymax) + 1/interval, 1/interval)
    xlocs = np.arange(np.floor(gxmin), np.ceil(gxmax) + 1/interval, 1/interval)
    if crosses:
        xlocs[xlocs > 180] -= 360

    gl.xlocator = mticker.FixedLocator(xlocs)
    gl.ylocator = mticker.FixedLocator(ylocs)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}


def get_shaded(ptopo, contour_colormap):
    maxvalue = contour_colormap.vmax
    ls1 = LightSource(azdeg=300, altdeg=45)
    ls2 = LightSource(azdeg=45, altdeg=45)
    intensity1 = ls1.hillshade(ptopo, fraction=0.25, vert_exag=VERT_EXAG)
    intensity2 = ls2.hillshade(ptopo, fraction=0.25, vert_exag=VERT_EXAG)
    intensity = intensity1 * 0.5 + intensity2 * 0.5

    ptoposc = ptopo / maxvalue
    rgba = contour_colormap.cmap(ptoposc)
    rgb = np.squeeze(rgba)

    draped_hsv = ls1.blend_hsv(rgb, np.expand_dims(intensity, 2))

    return draped_hsv


def draw_title(imt, container, operator):
    """Draw the map title.
    Args:
        imt (str): IMT that is being drawn on the map ('MMI', 'PGV',
            'PGA', 'SA(x.y)').
        isContour (bool): If true, use input imt, otherwise use MMI.
    """
    # Add a title
    edict = container.getMetadata()['input']['event_information']
    hlon = float(edict['longitude'])
    hlat = float(edict['latitude'])
    eloc = edict['event_description']
    try:
        etime = datetime.strptime(edict['origin_time'],
                                  constants.TIMEFMT)
    except ValueError:
        etime = datetime.strptime(edict['origin_time'],
                                  constants.ALT_TIMEFMT)
    timestr = etime.strftime('%b %d, %Y %H:%M:%S')
    mag = float(edict['magnitude'])
    if hlon < 0:
        lonstr = 'W%.2f' % np.abs(hlon)
    else:
        lonstr = 'E%.2f' % hlon
    if hlat < 0:
        latstr = 'S%.2f' % np.abs(hlat)
    else:
        latstr = 'N%.2f' % hlat
    dep = float(edict['depth'])
    eid = edict['event_id']
    if len(eid) <= 10:
        fmt = ('%s ShakeMap (%s): %s\n %s UTC M%.1f %s %s '
               'Depth: %.1fkm ID:%s')
    else:
        fmt = ('%s ShakeMap (%s): %s\n %s UTC M%.1f %s %s '
               'Depth: %.1fkm\nID:%s')
    tstr = fmt % (operator, imt, eloc, timestr, mag, latstr,
                  lonstr, dep, eid)
    plt.title(tstr, fontsize=10, verticalalignment='bottom')


def draw_stations(ax, stations, imt, intensity_colormap, geoproj, fill=True):
    """Draw station locations on the map.
    Args:
        ax (matplotlib Axes): Axes on which to draw stations.
        fill (bool): Whether or not to fill symbols.
        imt (str): One of ('MMI', 'PGA', 'PGV', or 'SA(x.x)')
    """
    dimt = imt.lower()

    # get the locations and values of the MMI observations
    mmi_dict = {
        'lat': [],
        'lon': [],
        'mmi': []
    }
    inst_dict = {
        'lat': [],
        'lon': [],
        'mmi': [],
        dimt: []
    }
    # Get the locations and values of the observed/instrumented
    # observations.
    for feature in stations['features']:
        lon, lat = feature['geometry']['coordinates']
        net = feature['properties']['network'].lower()

        # If the network matches one of these then it is an MMI
        # observation
        if net in ['dyfi', 'mmi', 'intensity', 'ciim']:
            # Append data from MMI features
            amplitude = feature['properties']['intensity']
            if amplitude == 'null':
                amplitude = 'nan'
            mmi_dict['mmi'].append(float(amplitude))
            mmi_dict['lat'].append(lat)
            mmi_dict['lon'].append(lon)
        else:
            # Otherwise, the feature is an instrument

            # If dimt is MMI then we have to use the converted value
            # from an instrumental value
            try:
                mmi_conv = float(feature['properties']['intensity'])
            except ValueError:
                mmi_conv = np.nan
            if dimt == 'mmi':
                inst_dict[dimt].append(mmi_conv)
                inst_dict['lat'].append(lat)
                inst_dict['lon'].append(lon)
            else:
                # Other IMTs are given in the channel for instruments
                channel = feature['properties']['channels'][0]
                for amplitude in channel['amplitudes']:
                    if amplitude['name'] != dimt:
                        continue
                    if amplitude['value'] == 'null':
                        continue
                    inst_dict[dimt].append(float(amplitude['value']))
                    inst_dict['lat'].append(lat)
                    inst_dict['lon'].append(lon)
                    # Add mmi also for the fill color of symbols
                    inst_dict['mmi'].append(mmi_conv)

    if fill:
        # Fill the symbols with the color for the intensity
        for i in range(len(mmi_dict['lat'])):
            mlat = mmi_dict['lat'][i]
            mlon = mmi_dict['lon'][i]
            mmi = mmi_dict['mmi'][i]
            mcolor = intensity_colormap.getDataColor(mmi)
            ax.plot(mlon, mlat, 'o', markerfacecolor=mcolor,
                    markeredgecolor='k', markersize=4, mew=0.5,
                    zorder=STATIONS_ZORDER, transform=geoproj)

        for i in range(len(inst_dict['lat'])):
            mlat = inst_dict['lat'][i]
            mlon = inst_dict['lon'][i]
            #
            # TODO: Make the fill color correspond to the mmi
            # obtained from the IMT.
            #
            mmi = inst_dict['mmi'][i]
            mcolor = intensity_colormap.getDataColor(mmi)
            ax.plot(mlon, mlat, '^',
                    markerfacecolor=mcolor, markeredgecolor='k',
                    markersize=6, zorder=STATIONS_ZORDER, mew=0.5,
                    transform=geoproj)
    else:
        # Do not fill symbols

        # plot MMI as small circles
        mmilat = mmi_dict['lat']
        mmilon = mmi_dict['lon']
        ax.plot(mmilon, mmilat, 'ko', fillstyle='none', mew=0.5,
                markersize=4, zorder=STATIONS_ZORDER, transform=geoproj)

        # plot MMI as slightly larger triangles
        instlat = inst_dict['lat']
        instlon = inst_dict['lon']
        ax.plot(instlon, instlat, 'k^', fillstyle='none', mew=0.5,
                markersize=6, zorder=STATIONS_ZORDER, transform=geoproj)


def get_draped(data, topodata, colormap):
    """Get array of data "draped" on topography.
    Args:
        data (ndarray): 2D Numpy array.
        topodata (ndarray): 2D Numpy array.
    Returns:
        ndarray: Numpy array of data draped on topography.
    """

    maxvalue = colormap.vmax
    mmisc = data / maxvalue
    rgba_img = colormap.cmap(mmisc)

    rgb = np.squeeze(rgba_img[:, :, 0:3])
    # use lightsource class to make our shaded topography
    ls = LightSource(azdeg=315, altdeg=45)

    ls1 = LightSource(azdeg=300, altdeg=45)
    ls2 = LightSource(azdeg=45, altdeg=45)
    intensity1 = ls1.hillshade(topodata, fraction=0.25,
                               vert_exag=VERT_EXAG)

    intensity2 = ls2.hillshade(topodata, fraction=0.25,
                               vert_exag=VERT_EXAG)

    intensity = intensity1 * 0.5 + intensity2 * 0.5

    draped_hsv = ls.blend_hsv(rgb, np.expand_dims(intensity, 2))

    return draped_hsv


def _clip_bounds(bbox, filename):
    """Clip input fiona-compatible vector file to input bounding box.

    :param bbox:
      Tuple of (xmin,ymin,xmax,ymax) desired clipping bounds.
    :param filename:
      Input name of file containing vector data in a format compatible
      with fiona.
    :returns:
      Shapely Geometry object (Polygon or MultiPolygon).
    """
    f = fiona.open(filename, 'r')
    shapes = list(f.items(bbox=bbox))
    xmin, ymin, xmax, ymax = bbox
    newshapes = []
    bboxpoly = sPolygon([(xmin, ymax), (xmax, ymax),
                         (xmax, ymin), (xmin, ymin), (xmin, ymax)])
    for tshape in shapes:
        myshape = sShape(tshape[1]['geometry'])
        intshape = myshape.intersection(bboxpoly)
        newshapes.append(intshape)
        newshapes.append(myshape)
    gc = GeometryCollection(newshapes)
    f.close()
    return gc


def draw_intensity(container, topobase, oceanfile, outpath, operator,
                   borderfile=None, override_scenario=False):
    """Create a contour map showing MMI contours over greyscale population.

    :param shakegrid:
      ShakeGrid object.
    :param popgrid:
      Grid2D object containing population data.
    :param oceanfile:
      String path to file containing ocean vector data in a format compatible
      with fiona.
    :param oceangridfile:
      String path to file containing ocean grid data .
    :param outpath:
      String path where output intensity.pdf and intensity.jpg files will
      be made.
    :param make_png:
      Boolean indicating whether a PNG version of the file should also be
      created in the same output folder as the PDF.
    :returns:
      Tuple containing:
        - Name of PNG file created, or None if PNG output not specified.
        - Cities object containing the cities that were rendered on the
          contour map.
    """
    # get the geodict for the ShakeMap
    comp = container.getComponents('MMI')[0]
    imtdict = container.getIMTGrids('MMI', comp)
    mmigrid = imtdict['mean']

    gd = mmigrid.getGeoDict()

    # Retrieve the epicenter - this will get used on the map
    origin = container.getRuptureObject().getOrigin()
    center_lat = origin.lat
    center_lon = origin.lon

    # load the cities data, limit to cities within shakemap bounds
    allcities = Cities.fromDefault()
    cities = allcities.limitByBounds((gd.xmin, gd.xmax, gd.ymin, gd.ymax))

    # get the map boundaries and figure size
    bounds, figsize = get_map_info(gd, center_lat)

    # Create the MercatorMap object, which holds a separate but identical
    # axes object used to determine collisions between city labels.
    mmap = MercatorMap(bounds, figsize, cities, padding=0.5,
                       dimensions=[0.1, 0.1, 0.8, 0.8])
    fig = mmap.figure
    ax = mmap.axes

    # Draw the map scale in the unoccupied lower corner.img_e
    corner = 'll'
    draw_scale(ax, corner, pady=0.05, padx=0.05, zorder=SCALE_ZORDER)

    # this needs to be done here so that city label collision
    # detection will work
    fig.canvas.draw()

    # get the geographic projection object
    geoproj = mmap.geoproj
    # get the mercator projection object
    proj = mmap.proj
    # get the proj4 string - used by Grid2D project() method
    projstr = proj.proj4_init

    # get the projected MMI and topo grids
    pmmigrid, ptopogrid = get_projected_grids(mmigrid, topobase, projstr)

    # get the projected geodict
    proj_gd = pmmigrid.getGeoDict()

    # Use our GMT-inspired palette class to create population and MMI colormaps
    mmimap = ColorPalette.fromPreset('mmi')

    # drape the intensity data over the topography
    pmmi_data = pmmigrid.getData()
    ptopo_data = ptopogrid.getData()
    draped_hsv = get_draped(pmmi_data, ptopo_data, mmimap)

    # set the image extent to that of the data (do we need this?)
    # set the axes border width
    plt.sca(ax)
    ax.set_xlim(proj_gd.xmin, proj_gd.xmax)
    ax.set_ylim(proj_gd.ymin, proj_gd.ymax)
    img_extent = (proj_gd.xmin, proj_gd.xmax, proj_gd.ymin, proj_gd.ymax)
    plt.imshow(draped_hsv, origin='upper', extent=img_extent,
               zorder=IMG_ZORDER, interpolation='none')

    # draw 10m res coastlines
    ax.coastlines(resolution="10m", zorder=BORDER_ZORDER)

    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    ax.add_feature(states_provinces, edgecolor='black', zorder=BORDER_ZORDER)

    # draw country borders using natural earth data set
    if borderfile is not None:
        borders = ShapelyFeature(Reader(borderfile).geometries(),
                                 ccrs.PlateCarree())
        ax.add_feature(borders, zorder=BORDER_ZORDER,
                       edgecolor='black', linewidth=2, facecolor='none')

    # clip the ocean data to the shakemap
    bbox = (gd.xmin, gd.ymin, gd.xmax, gd.ymax)
    oceanshapes = _clip_bounds(bbox, oceanfile)

    ax.add_feature(ShapelyFeature(oceanshapes, crs=geoproj),
                   facecolor=WATERCOLOR, zorder=OCEAN_ZORDER)

    # draw cities
    mmap.drawCities(shadow=True, zorder=CITIES_ZORDER,
                    draw_dots=True)

    # draw the station data on the map
    stations = container.getStationDict()

    draw_stations(ax, stations, 'MMI', mmimap, geoproj)

    # Draw the epicenter as a black star
    plt.sca(ax)
    plt.plot(center_lon, center_lat, 'k*', markersize=16,
             zorder=EPICENTER_ZORDER, transform=geoproj)

    # draw the map title
    draw_title('MMI', container, operator)

    # is this event a scenario?
    is_scenario = True

    if is_scenario:
        plt.text(center_lon, center_lat, 'SCENARIO', fontsize=64,
                 zorder=SCENARIO_ZORDER, transform=geoproj,
                 alpha=0.2, color='red', horizontalalignment='center')

    # draw the rupture polygon(s) in black, if not point rupture
    rupture = container.getRuptureObject()
    if not isinstance(rupture, PointRupture):
        lats = rupture.lats
        lons = rupture.lons
        ax.plot(lons, lats, 'k', lw=2, zorder=FAULT_ZORDER, transform=geoproj)

    # draw graticules, ticks, tick labels
    draw_graticules(ax, *bounds)

    # draw a separate intensity colorbar in a separate axes
    draw_colorbar(fig, mmimap)

    # # make the map border thicker
    # the left/top edges don't line up here
    # plt.sca(ax)
    # lw = 5.0
    # ax.outline_patch.set_zorder(AXES_ZORDER)
    # ax.outline_patch.set_linewidth(lw)

    # create pdf and png output file names
    pdf_file = os.path.join(outpath, 'intensity.pdf')
    jpg_file = os.path.join(outpath, 'intensity.jpg')

    # save to pdf/jpeg
    plt.savefig(pdf_file)
    plt.savefig(jpg_file)

    return (pdf_file, jpg_file)


def draw_contour(container, imtype, topobase, oceanfile, outpath,
                 operator, filter_size, borderfile=None, is_scenario=False):

    comp = container.getComponents(imtype)[0]
    imtdict = container.getIMTGrids(imtype, comp)
    imtgrid = imtdict['mean']

    gd = imtgrid.getGeoDict()

    # Retrieve the epicenter - this will get used on the map
    origin = container.getRuptureObject().getOrigin()
    center_lat = origin.lat

    # load the cities data, limit to cities within shakemap bounds
    allcities = Cities.fromDefault()
    cities = allcities.limitByBounds((gd.xmin, gd.xmax, gd.ymin, gd.ymax))

    # get the map boundaries and figure size
    bounds, figsize = get_map_info(gd, center_lat)

    # Create the MercatorMap object, which holds a separate but identical
    # axes object used to determine collisions between city labels.
    mmap = MercatorMap(bounds, figsize, cities, padding=0.5,
                       dimensions=[0.1, 0.1, 0.8, 0.8])
    fig = mmap.figure
    ax = mmap.axes
    # this needs to be done here so that city label collision
    # detection will work
    fig.canvas.draw()

    # get the geographic projection object
    geoproj = mmap.geoproj
    # get the mercator projection object
    proj = mmap.proj
    # get the proj4 string - used by Grid2D project() method
    projstr = proj.proj4_init

    # get the projected MMI and topo grids
    pimtgrid, ptopogrid = get_projected_grids(imtgrid, topobase, projstr)

    # get the projected geodict
    proj_gd = pimtgrid.getGeoDict()

    # convert units if necessary
    pimtdata = pimtgrid.getData()
    if imtype == 'MMI':
        pass
    elif imtype == 'PGV':
        pimtdata = np.exp(pimtdata)
    else:
        pimtdata = np.exp(pimtdata) * 100

    # get the draped topo data
    topo_colormap = ColorPalette.fromPreset('shaketopo')
    hillshade = get_shaded(ptopogrid.getData(), topo_colormap)

    # draw the draped intensity data
    plt.sca(ax)
    ax.set_xlim(proj_gd.xmin, proj_gd.xmax)
    ax.set_ylim(proj_gd.ymin, proj_gd.ymax)
    img_extent = (proj_gd.xmin, proj_gd.xmax, proj_gd.ymin, proj_gd.ymax)
    plt.imshow(hillshade, origin='upper', extent=img_extent,
               zorder=IMG_ZORDER, interpolation='none')

    config = container.getConfig()
    gmice = get_object_from_config('gmice', 'modeling', config)
    gmice_imts = gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES
    gmice_pers = gmice.DEFINED_FOR_SA_PERIODS

    oqimt = imt.from_string(imtype)
    component = container.getComponents(imtype)[0]

    if imtype == 'MMI' or not isinstance(oqimt, tuple(gmice_imts)) or \
            (isinstance(oqimt, imt.SA) and oqimt.period not in gmice_pers):
        my_gmice = None
    else:
        my_gmice = gmice

    # call the contour module in plotting to get the vertices of the contour
    # lines
    contour_objects = contour(container, imtype,
                              component, filter_size,
                              my_gmice)

    # cartopy shapely feature has some weird behaviors, so I had to go rogue
    # and draw contour lines/labels myself.
    # draw dashed contours first, the ones over land will be overridden by
    # solid contours
    npoints = []
    for contour_object in contour_objects:
        props = contour_object['properties']
        lw = props['weight']
        multi_lines = sShape(contour_object['geometry'])
        pmulti_lines = proj.project_geometry(multi_lines, src_crs=geoproj)
        for multi_line in pmulti_lines:
            pmulti_line = mapping(multi_line)['coordinates']
            x, y = zip(*pmulti_line)
            npoints.append(len(x))
            ax.plot(x, y, color=props['color'], linestyle='dashed',
                    zorder=DASHED_CONTOUR_ZORDER)

    white_box = dict(boxstyle="round",
                     ec=(0, 0, 0),
                     fc=(1., 1, 1),
                     color='k'
                     )

    # only label lines with lots of points
    npoints = np.array(npoints)
    # min_npoints = npoints.mean() - (npoints.std()/2)
    min_npoints = npoints.mean()

    # draw solid contours next - the ones over water will be covered by ocean
    # polygon
    for contour_object in contour_objects:
        props = contour_object['properties']
        multi_lines = sShape(contour_object['geometry'])
        pmulti_lines = proj.project_geometry(multi_lines, src_crs=geoproj)
        for multi_line in pmulti_lines:
            pmulti_line = mapping(multi_line)['coordinates']
            x, y = zip(*pmulti_line)
            ax.plot(x, y, color=props['color'], linestyle='solid',
                    zorder=CONTOUR_ZORDER)
            if len(x) > min_npoints:
                # try to label each segment with black text in a white box
                # xc, yc = choose_label_vertex(x, y, proj_gd.xmin, proj_gd.xmax,
                #                              proj_gd.ymin, proj_gd.ymax)
                xc = x[int(len(x)/3)]
                yc = y[int(len(y)/3)]
                if label_close_to_edge(xc, yc, proj_gd.xmin, proj_gd.xmax,
                                       proj_gd.ymin, proj_gd.ymax):
                    continue
                # TODO: figure out if box is going to go outside the map, if so
                # choose a different point on the line.
                ax.text(xc, yc, '%.1f' % props['value'], size=8,
                        ha="center", va="center",
                        bbox=white_box, zorder=AXES_ZORDER-1)

    ax.outline_patch.set_zorder(AXES_ZORDER)

    # clip the ocean data to the shakemap
    bbox = (gd.xmin, gd.ymin, gd.xmax, gd.ymax)
    oceanshapes = _clip_bounds(bbox, oceanfile)

    ax.add_feature(ShapelyFeature(oceanshapes, crs=geoproj),
                   facecolor=WATERCOLOR, zorder=OCEAN_ZORDER)

    # draw 10m res coastlines
    ax.coastlines(resolution="10m", zorder=COAST_ZORDER, linewidth=3)

    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    ax.add_feature(states_provinces, edgecolor='black', zorder=COAST_ZORDER)

    # draw graticules, ticks, tick labels
    draw_graticules(ax, *bounds)

    # Draw the map scale in the unoccupied lower corner.
    corner = 'll'
    draw_scale(ax, corner, pady=0.05, padx=0.05, zorder=SCALE_ZORDER)

    # draw cities
    mmap.drawCities(shadow=True, zorder=CITIES_ZORDER,
                    draw_dots=True)

    # draw the rupture polygon(s) in black, if not point rupture
    rupture = container.getRuptureObject()
    if not isinstance(rupture, PointRupture):
        lats = rupture.lats
        lons = rupture.lons
        ax.plot(lons, lats, 'k', lw=2, zorder=FAULT_ZORDER, transform=geoproj)

    # draw the station data on the map
    stations = container.getStationDict()
    mmimap = ColorPalette.fromPreset('mmi')
    draw_stations(ax, stations, imtype, mmimap, geoproj)

    draw_title(imtype, container, operator)

    # draw a separate intensity colorbar in a separate axes
    draw_colorbar(fig, mmimap)

    # save plot to file
    fileimt = oq_to_file(imtype)
    plt.draw()
    pdf_file = os.path.join(outpath, '%s_contour.pdf' % (fileimt))
    jpg_file = os.path.join(outpath, '%s_contour.jpg' % (fileimt))
    plt.savefig(pdf_file)
    plt.savefig(jpg_file)

    return (pdf_file, jpg_file)
