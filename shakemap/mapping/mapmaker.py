# stdlib imports
from datetime import datetime

# third party imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as ccrs  # projections
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import ShapelyFeature
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader

from shapely.geometry import shape as sShape
from shapely.geometry import Polygon as sPolygon
from shapely.geometry import GeometryCollection

import pyproj
import fiona

# neic imports
from mapio.gmt import GMTGrid
from impactutils.mapping.mercatormap import MercatorMap
from impactutils.mapping.city import Cities
from impactutils.colors.cpalette import ColorPalette
from impactutils.mapping.scalebar import draw_scale

# local imports
from shakelib.rupture.point_rupture import PointRupture
from shakelib.rupture import constants

# define some constants
WATERCOLOR = '#7AA1DA'
FIGWIDTH = 7.0
FILTER_SMOOTH = 5.0
XOFFSET = 4  # how many pixels between the city dot and the city text
VERT_EXAG = 0.001  # what is the vertical exaggeration for hillshade

# define the zorder values for various map components
# all of the zorder values for different plotted parameters
# elements with a higher zorder will plot on top of lower elements.
IMG_ZORDER = 1
ROAD_ZORDER = 5
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

# default font for cities
# DEFAULT_FONT = 'DejaVu Sans'
DEFAULT_FONT = 'Bitstream Vera Sans'

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
    fmt = ('%s ShakeMap (%s): %s\n %s UTC M%.1f %s %s '
           'Depth: %.1fkm ID:%s')
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
            channel = feature['properties']['channels'][0]
            for amplitude in channel['amplitudes']:
                if amplitude['name'] != 'mmi':
                    continue
                mmi_dict['mmi'].append(float(amplitude['value']))
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
    ls = LightSource(azdeg=135, altdeg=45)

    ls1 = LightSource(azdeg=120, altdeg=45)
    ls2 = LightSource(azdeg=225, altdeg=45)
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


def draw_intensity(container, topofile, oceanfile, basename, operator,
                   borderfile=None, is_scenario=False):
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
    :param basename:
      String path containing desired output PDF base name, i.e.,
      /home/pager/exposure.  ".pdf" and ".png" files will
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
    smdict = mmigrid.getGeoDict()

    gd = mmigrid.getGeoDict()

    # Retrieve the epicenter - this will get used on the map
    origin = container.getRuptureObject().getOrigin()
    center_lat = origin.lat
    center_lon = origin.lon

    # load the cities data, limit to cities within shakemap bounds
    allcities = Cities.fromDefault()
    cities = allcities.limitByBounds((gd.xmin, gd.xmax, gd.ymin, gd.ymax))

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
    bbox = (xmin, ymin, xmax, ymax)
    bounds = (xmin, xmax, ymin, ymax)
    figsize = (FIGWIDTH, figheight)

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

    # resample shakemap to topogrid
    # get the geodict for the topo file
    topodict = GMTGrid.getFileGeoDict(topofile)[0]

    # get a geodict that is aligned with topo, but inside shakemap
    sampledict = topodict.getBoundsWithin(smdict)

    mmigrid = mmigrid.interpolateToGrid(sampledict)

    gd = mmigrid.getGeoDict()

    # get topo layer and project it
    topogrid = GMTGrid.load(topofile, samplegeodict=sampledict, resample=False)
    ptopogrid = topogrid.project(projstr)

    # project the MMI data
    pmmigrid = mmigrid.project(projstr)

    # get the projected geodict
    proj_gd = pmmigrid.getGeoDict()

    # Use our GMT-inspired palette class to create population and MMI colormaps
    mmimap = ColorPalette.fromPreset('mmi')

    # drape the intensity data over the topography
    pmmi_data = pmmigrid.getData()
    ptopo_data = ptopogrid.getData()
    draped_hsv = get_draped(pmmi_data, ptopo_data, mmimap)

    # set the image extent to that of the data (do we need this?)
    img_extent = (proj_gd.xmin, proj_gd.xmax, proj_gd.ymin, proj_gd.ymax)
    draped_img = plt.imshow(draped_hsv, origin='upper', extent=img_extent,
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

    # draw meridians and parallels using Cartopy's functions for that
    gl = ax.gridlines(draw_labels=True,
                      linewidth=2, color=(0.9, 0.9, 0.9),
                      alpha=0.5, linestyle='-',
                      zorder=GRATICULE_ZORDER)
    gl.xlabels_top = False
    gl.xlabels_bottom = False
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.xlines = True

    # let's floor/ceil the edges to nearest half a degree
    gxmin = np.floor(xmin * 2) / 2
    gxmax = np.ceil(xmax * 2) / 2
    gymin = np.floor(ymin * 2) / 2
    gymax = np.ceil(ymax * 2) / 2

    # shakemap way
    ylocs = np.arange(np.ceil(gymin), np.floor(gymax) + 1, 1.0)
    xlocs = np.arange(np.ceil(gxmin), np.floor(gxmax) + 1, 1.0)

    # PAGER way
    # xlocs = np.linspace(gxmin, gxmax+0.5, num=5)
    # ylocs = np.linspace(gymin, gymax+0.5, num=5)

    gl.xlocator = mticker.FixedLocator(xlocs)
    gl.ylocator = mticker.FixedLocator(ylocs)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'black'}
    gl.ylabel_style = {'size': 15, 'color': 'black'}

    # TODO - figure out x/y axes data coordinates
    # corresponding to 10% from left and 10% from top
    # use geoproj and proj
    dleft = 0.01
    dtop = 0.97
    proj_str = proj.proj4_init
    merc_to_dd = pyproj.Proj(proj_str)

    # use built-in transforms to get from axes units to data units
    display_to_data = ax.transData.inverted()
    axes_to_display = ax.transAxes

    # these are x,y coordinates in projected space
    yleft, t1 = display_to_data.transform(
        axes_to_display.transform((dleft, 0.5)))
    t2, xtop = display_to_data.transform(
        axes_to_display.transform((0.5, dtop)))

    # these are coordinates in lon,lat space
    yleft_dd, t1_dd = merc_to_dd(yleft, t1, inverse=True)
    t2_dd, xtop_dd = merc_to_dd(t2, xtop, inverse=True)

    # drawing our own tick labels INSIDE the plot, as
    # Cartopy doesn't seem to support this.
    yrange = ymax - ymin
    xrange = xmax - xmin
    ddlabelsize = 12
    for xloc in gl.xlocator.locs:
        outside = xloc < xmin or xloc > xmax
        # don't draw labels when we're too close to either edge
        near_edge = (xloc-xmin) < (xrange*0.1) or (xmax-xloc) < (xrange*0.1)
        if outside or near_edge:
            continue
        xtext = r'$%.1f^\circ$W' % (abs(xloc))
        ax.text(xloc, xtop_dd, xtext,
                fontsize=ddlabelsize, zorder=GRATICULE_ZORDER, ha='center',
                fontname=DEFAULT_FONT,
                transform=ccrs.Geodetic())

    for yloc in gl.ylocator.locs:
        outside = yloc < gd.ymin or yloc > gd.ymax
        # don't draw labels when we're too close to either edge
        near_edge = (yloc-gd.ymin) < (yrange *
                                      0.1) or (gd.ymax-yloc) < (yrange*0.1)
        if outside or near_edge:
            continue
        if yloc < 0:
            ytext = r'$%.1f^\circ$S' % (abs(yloc))
        else:
            ytext = r'$%.1f^\circ$N' % (abs(yloc))
        ax.text(yleft_dd, yloc, ytext,
                fontsize=ddlabelsize, zorder=GRATICULE_ZORDER, va='center',
                fontname=DEFAULT_FONT,
                transform=ccrs.Geodetic())

    # draw cities
    mapcities = mmap.drawCities(shadow=False, zorder=CITIES_ZORDER)

    # draw the station data on the map
    stations = container.getStationDict()

    draw_stations(ax, stations, 'MMI', mmimap, geoproj)

    # Draw the map scale in the unoccupied lower corner.
    corner = 'll'
    draw_scale(ax, corner, pady=0.05, padx=0.05, zorder=SCALE_ZORDER)

    # Draw the epicenter as a black star
    plt.sca(ax)
    plt.plot(center_lon, center_lat, 'k*', markersize=16,
             zorder=EPICENTER_ZORDER, transform=geoproj)

    # draw the map title
    draw_title('MMI', container, operator)

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

    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("bottom", size="5%",
    #                           pad=0.05, map_projection=proj)
    cax = fig.add_axes([0.1, 0.035, 0.8, 0.050])

    # making up our own colorbar object here because the default
    # pyplot functionality doesn't seem to do the job.

    # sm = plt.cm.ScalarMappable(cmap=mmimap.cmap,
    #                            norm=plt.Normalize(vmin=1, vmax=10))
    # sm._A = []
    # plt.colorbar(sm, orientation='horizontal',
    #              shrink=0.8, fraction=0.05, cax=cax)
    # plt.colorbar(draped_img, orientation='horizontal',
    #              shrink=0.80, fraction=0.05, cax=cax)

    # create pdf and png output file names
    pdf_file = basename+'.pdf'
    png_file = basename+'.png'

    # save to pdf
    plt.savefig(pdf_file)
    plt.savefig(png_file)

    return (pdf_file, png_file, mapcities)
