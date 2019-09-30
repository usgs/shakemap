# stdlib imports
from datetime import datetime
from collections import defaultdict

# third party imports
import numpy as np
import matplotlib.image as image
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LightSource
import matplotlib.patheffects as path_effects
from matplotlib.font_manager import FontProperties
from matplotlib import patches

import cartopy.crs as ccrs  # projections
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pyproj

from shapely.geometry import shape as sShape
from shapely.geometry import Polygon as sPolygon
from shapely.geometry import LineString as sLineString
from shapely.geometry import GeometryCollection
from shapely.geometry import mapping

import fiona
from openquake.hazardlib import imt

# neic imports
from impactutils.mapping.mercatormap import MercatorMap
from impactutils.colors.cpalette import ColorPalette
from impactutils.mapping.scalebar import draw_scale
from impactutils.textformat.text import set_num_precision
from mapio.grid2d import Grid2D
from mapio.geodict import GeoDict


# local imports
from shakelib.rupture.point_rupture import PointRupture
from shakelib.rupture import constants
from shakelib.rupture.factory import rupture_from_dict
from shakelib.plotting.contour import contour, getContourLevels
from shakelib.gmice.wgrw12 import WGRW12
from shakemap.utils.utils import get_object_from_config

# define some constants
WATERCOLOR = '#7AA1DA'
FIGWIDTH = 9.5
FIGHEIGHT = 10.0
XOFFSET = 4  # how many pixels between the city dot and the city text
VERT_EXAG = 0.1  # what is the vertical exaggeration for hillshade

# define the zorder values for various map components
# all of the zorder values for different plotted parameters
# elements with a higher zorder will plot on top of lower elements.
IMG_ZORDER = 1
ROAD_ZORDER = 5000
COAST_ZORDER = 11
CONTOUR_ZORDER = 800
DASHED_CONTOUR_ZORDER = 1002
OCEAN_ZORDER = 1000
EPICENTER_ZORDER = 1100
BORDER_ZORDER = 1110
STATIONS_ZORDER = 1150
FAULT_ZORDER = 1160
CITIES_ZORDER = 1200
GRATICULE_ZORDER = 1200
SCALE_ZORDER = 1500
SCENARIO_ZORDER = 2000
AXES_ZORDER = 3000

# default font for cities
# DEFAULT_FONT = 'DejaVu Sans'
DEFAULT_FONT = 'Bitstream Vera Sans'

# what color should the scenario watermark be?
WATERMARK_COLOR = (0.85, 0.85, 0.85)
WATERMARK_ALPHA = 0.4

# what fraction of the width/height of the map can the contour labels be
# from the edge before we decide we don't want to draw it there?
MAP_FRAC = 10

# define dictionary of MMI integer values to Roman numeral equivalents
MMI_LABELS = {
    '1': 'I',
    '2': 'II',
    '3': 'III',
    '4': 'IV',
    '5': 'V',
    '6': 'VI',
    '7': 'VII',
    '8': 'VIII',
    '9': 'IX',
    '10': 'X'
}

IMT_RANGES = {
    'PGV': (1e-2, 500),
    'PGA': (1e-4, 500),
    'SA(0.3)': (1e-4, 500),
    'SA(1.0)': (1e-4, 500),
    'SA(3.0)': (1e-4, 400)
}

DEG2KM = 111.191

# what is the colormap we want to use for non-MMI intensity values?
IMT_CMAP = 'PuRd'


def to_precision(x, p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(np.log10(x))
    tens = np.power(10, e - p + 1)
    n = np.floor(x/tens)

    if n < np.power(10, p - 1):
        e = e - 1
        tens = np.power(10, e - p+1)
        n = np.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens - x):
        n = n + 1

    if n >= np.power(10, p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p - 1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)


def _create_palette(imtype, levels):
    """Create a ColorPalette object from given levels and IMT type.

    Args:
        imtype (str): One of 'PGV','PGA','SA(0.3)',etc.
        levels (sequence): Sequence of contour levels.
    Returns:
        ColorPalette: ColorPalette using range of input data and IMT_CMAP.
    """
    # this method assumes that levels are in logspace
    if len(levels) > 1:
        if len(levels) % 2:
            levels.append(levels[-1])
        nsteps = 256
        z0 = np.linspace(np.log(levels[0]), np.log(levels[-2]), nsteps)
        z1 = np.linspace(np.log(levels[1]), np.log(levels[-1]), nsteps)
    else:
        z0 = np.array([levels[0], levels[0]*10])
        z1 = np.array([levels[0], levels[0]*10])
    cmap = plt.get_cmap(IMT_CMAP)
    palette = ColorPalette.fromColorMap(imtype, z0, z1, cmap, is_log=True)
    return palette


def _get_projected_grids(imtgrid, topobase, projstr):
    """Resample IMT grid to topography, return projected versions of each.

    Args:
        imtgrid (Grid2D): IMT mean grid.
        topobase (Grid2D): Topography grid to trim and project.
        projstr (str): Proj4 initialization string.
    Returns:
        Grid2D: projected IMT mean grid, resampled to topography.
        Grid2D: projected topography grid, trimmed to size of IMT.
    """
    # resample shakemap to topogrid
    # get the geodict for the topo file
    topodict = topobase.getGeoDict()

    # get a geodict that is aligned with topo, but inside shakemap
    sampledict = topodict.getBoundsWithin(imtgrid.getGeoDict())

    simtgrid = imtgrid.interpolateToGrid(sampledict, method='nearest')

    # get topo layer and project it
    topogrid = topobase.interpolateToGrid(sampledict, method='linear')

    # resampling 32bit floats gives odd results... upcasting to 64bit
    topogrid._data = topogrid._data.astype(np.float64)

    # project the topography data
    ptopogrid = topogrid.project(projstr, method='bilinear')

    # project the MMI data
    pimtgrid = simtgrid.project(projstr)

    return (pimtgrid, ptopogrid)


def _get_map_info(gd):
    """Get the desired bounds of the map and the figure size (in).

    Args:
        gd (GeoDict): GeoDict to use as basis for map boundaries.

    Returns:
        tuple: xmin,xmax,ymin,ymax Extent of map.
        tuple: Figure size (width,height) in inches.
    """
    # define the map
    # first cope with stupid 180 meridian
    xmin, xmax, ymin, ymax = (gd.xmin, gd.xmax, gd.ymin, gd.ymax)
    if xmin > xmax:
        xmax = xmax + 360

    center_lon = (xmin + xmax)/2.0

    proj = ccrs.Mercator(central_longitude=center_lon,
                         min_latitude=ymin,
                         max_latitude=ymax,
                         globe=None)
    pproj = pyproj.Proj(proj.proj4_init)
    pxmin, pymin = pproj(xmin, ymin)
    pxmax, pymax = pproj(xmax, ymax)
    pwidth = pxmax - pxmin
    pheight = pymax - pymin

    bounds = (xmin, xmax, ymin, ymax)

    # Map aspect
    aspect = pwidth/pheight
    # This all seems unnecessary
    # fig_aspect = 1.0/(0.19 + 0.8/aspect)
    # figheight = FIGWIDTH / fig_aspect
    # figsize = (FIGWIDTH, figheight)
    # Make the figsize a constant, other functions will fit the map and
    # legends within the available space
    figsize = (FIGWIDTH, FIGHEIGHT)
    return (bounds, figsize, aspect)


def _draw_imt_legend(fig, palette, imtype, gmice, process_time, map_version,
                     point_source, tdict):
    """Create a legend axis for non MMI plots.

    Args:
        fig (Figure): Matplotlib Figure object.
        levels (sequence): Sequence of contour levels.
        palette (ColorPalette): ColorPalette using range of input data and
            IMT_CMAP.
        imtype (str): One of 'PGV','PGA','SA(0.3)',etc.
        gmice (GMICE object): The GMICE used for this map.
        process_time (str): The processing time of this map.
        map_version (str): The version of this map.
        point_source (bool): Is the rupture a point source?
        tdict (dict): Dictionary containing the text strings for printing
            on the maps (in the language of the user's choice).
    """
    imtlabel = imtype + ' ' + tdict['units'][imtype]
    # imtlabel = imtype

    cax = fig.add_axes([0.1, 0.13, 0.8, 0.02])
    plt.axis('off')
    cax_xmin, cax_xmax = cax.get_xlim()
    bottom, top = cax.get_ylim()
    plt.xlim(cax_xmin, cax_xmax)
    plt.ylim(bottom, top)

    firstcol_width = 0.15

    font0 = FontProperties()
    alignment = {
        'horizontalalignment': 'center',
        'verticalalignment': 'center'
    }
    font0.set_weight('bold')

    xloc = firstcol_width/2
    plt.text(xloc, 0.5, imtlabel,
             fontproperties=font0, **alignment)
    # draw top/bottom edges of table
    plt.plot([bottom, top], [bottom, bottom], 'k', clip_on=False)
    plt.plot([bottom, top], [top, top], 'k', clip_on=False)
    # draw left edge of table
    plt.plot([bottom, bottom], [bottom, top], 'k', clip_on=False)
    # draw right edge of first column
    plt.plot([firstcol_width, firstcol_width], [0, 1], 'k', clip_on=False)
    # draw right edge of table
    plt.plot([1, 1], [0, 1], 'k', clip_on=False)

    # get the MMI/IMT values we need
    itype = 'log'
    divisor = 1
    if imtype != 'PGV':
        divisor = 100
    dmin, dmax = IMT_RANGES[imtype]
    imt_values = np.log(getContourLevels(dmin, dmax, itype=itype)/divisor)
    if gmice.supports(imtype):
        mmi_values, _ = gmice.getMIfromGM(imt_values, imt.from_string(imtype))
    else:
        gmice = WGRW12()
        mmi_values, _ = gmice.getMIfromGM(imt_values, imt.from_string(imtype))
    mmi_colors = [palette.getDataColor(
        mmi, color_format='hex') for mmi in mmi_values]
    new_imts = []
    new_mmi_colors = []
    for mmic, imtv in zip(mmi_colors, imt_values):
        if mmic not in new_mmi_colors:
            new_imts.append(imtv)
            new_mmi_colors.append(mmic)

    width = (1 - firstcol_width)/len(new_imts)
    left = firstcol_width
    for mmic, imtv in zip(new_mmi_colors, imt_values):
        right = left + width
        px = [left, right, right, left, left]
        py = [top, top, bottom, bottom, top]
        plt.plot([right, right], [bottom, top], 'k')
        plt.fill(px, py, mmic, ec=mmic)
        xloc = left + width/2.0
        imtstr = "{0:.3g}".format(np.exp(imtv)*divisor)
        th = plt.text(xloc, 0.5, imtstr, fontproperties=font0, **alignment)
        th.set_path_effects(
            [path_effects.Stroke(linewidth=2.0,
                                 foreground='white'),
             path_effects.Normal()]
        )
        left = right

    # Explanation of symbols: triangle is instrument, circle is mmi,
    # epicenter is black star
    # thick black line is rupture (if available)
    cax = fig.add_axes([0.1, 0.09, 0.8, 0.04])
    plt.axis('off')
    cax_xmin, cax_xmax = cax.get_xlim()
    bottom, top = cax.get_ylim()
    plt.xlim(cax_xmin, cax_xmax)
    plt.ylim(bottom, top)
    item_sep = [0.2, 0.28, 0.15]
    left_offset = 0.005
    label_pad = 0.02

    yloc_sixth_row = 0.6
    yloc_seventh_row = 0.15

    # Instrument
    triangle_marker_x = left_offset
    triangle_text_x = triangle_marker_x + label_pad
    plt.plot(triangle_marker_x, yloc_seventh_row, '^', markerfacecolor='w',
             markeredgecolor='k', markersize=6, mew=0.5, clip_on=False)
    plt.text(triangle_text_x,
             yloc_seventh_row,
             tdict['legend']['instrument'],
             va='center',
             ha='left')

    # Macroseismic
    circle_marker_x = triangle_text_x + item_sep[0]
    circle_text_x = circle_marker_x + label_pad
    plt.plot(circle_marker_x,
             yloc_seventh_row, 'o',
             markerfacecolor='w',
             markeredgecolor='k',
             markersize=4,
             mew=0.5)
    plt.text(circle_text_x,
             yloc_seventh_row,
             tdict['legend']['intensity'],
             va='center',
             ha='left')

    # Epicenter
    star_marker_x = circle_marker_x + item_sep[1]
    star_text_x = star_marker_x + label_pad
    plt.plot(star_marker_x,
             yloc_seventh_row, 'k*',
             markersize=12,
             mew=0.5)
    plt.text(star_text_x,
             yloc_seventh_row,
             tdict['legend']['epicenter'],
             va='center',
             ha='left')

    if not point_source:
        rup_marker_x = star_marker_x + item_sep[2]
        rup_text_x = rup_marker_x + label_pad
        rwidth = 0.02
        rheight = 0.05
        rup = patches.Rectangle(
            xy=(rup_marker_x - rwidth,
                yloc_seventh_row-0.5*rheight),
            width=rwidth,
            height=rheight,
            linewidth=2,
            edgecolor='k',
            facecolor='w'
        )
        cax.add_patch(rup)
        plt.text(rup_text_x,
                 yloc_seventh_row,
                 tdict['legend']['rupture'],
                 va='center',
                 ha='left')

    # Add conversion reference and shakemap version/process time
    version_x = 1.0
    tpl = (tdict['legend']['version'], map_version,
           tdict['legend']['processed'], process_time)
    plt.text(version_x, yloc_sixth_row,
             '%s %i: %s %s' % tpl,
             ha='right', va='center')

    ref = gmice.name
    refx = 0
    plt.text(refx, yloc_sixth_row,
             '%s %s' % (tdict['legend']['scale'], ref),
             va='center')


def _draw_mmi_legend(fig, palette, gmice, process_time, map_version,
                     point_source, tdict):
    """Create a legend axis for MMI plots.

    Args:
        fig (Figure): Matplotlib Figure object.
        palette (ColorPalette): ColorPalette using range of input data and
            IMT_CMAP.
        gmice: A gmice object.
        process_time (str): Process time.
        map_version (int): ShakeMap version.
        point_source (bool): Is the rupture a PointRupture?
        tdict (dict): Dictionary containing the text strings for printing
            on the maps (in the language of the user's choice).

    """
    cax = fig.add_axes([0.1, 0.00, 0.8, 0.15])
    plt.axis('off')
    cax_xmin, cax_xmax = cax.get_xlim()
    bottom, top = cax.get_ylim()
    plt.xlim(cax_xmin, cax_xmax)
    plt.ylim(bottom, top)

    acceleration = [tdict['mmi_scale']['acc_label']]
    velocity = [tdict['mmi_scale']['vel_label']]

    imt_edges = np.array([0.5, 1.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    mmi_centers = np.array([1.0, 2.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    pga_values, _ = gmice.getGMfromMI(mmi_centers, imt.from_string('PGA'))
    pgv_values, _ = gmice.getGMfromMI(mmi_centers, imt.from_string('PGV'))
    pga_values = np.exp(pga_values)*100
    pgv_values = np.exp(pgv_values)
    pga_labels = ["{0:.3g}".format(set_num_precision(
        pga, 3, mode='float')) for pga in pga_values]
    pgv_labels = ["{0:.3g}".format(set_num_precision(
        pgv, 3, mode='float')) for pgv in pgv_values]

    pga_labels[0] = '<'+pga_labels[0]
    pga_labels[-1] = '>'+pga_labels[-1]
    pgv_labels[0] = '<'+pgv_labels[0]
    pgv_labels[-1] = '>'+pgv_labels[-1]

    acceleration += pga_labels
    velocity += pgv_labels

    yloc_first_row = 13/14
    yloc_second_row = 11/14
    yloc_third_row = 9/14
    yloc_fourth_row = 7/14
    yloc_fifth_row = 5/14
    yloc_sixth_row = 3/14
    yloc_seventh_row = 1.5/14

    yloc_first_line = 12/14
    yloc_second_line = 10/14
    yloc_third_line = 8/14
    yloc_fourth_line = 6/14
    # yloc_fifth_line = 4/14

    bottom = 4/14

    font0 = FontProperties()
    alignment = {
        'horizontalalignment': 'center',
        'verticalalignment': 'center'
    }
    font0.set_weight('bold')

    font1 = FontProperties()
    font1.set_weight('normal')

    # draw vertical cell separators
    sumwidth = 0.0
    gridleft = 0.0
    plt.plot([gridleft, gridleft], [bottom, top],
             'k', clip_on=False)  # left edge
    plt.plot([0, 1], [top, top], 'k', clip_on=False)
    plt.plot([0, 1], [bottom, bottom], 'k', clip_on=False)

    plt.plot([0, 1], [yloc_first_line, yloc_first_line],
             'k', clip_on=False)
    plt.plot([0, 1], [yloc_second_line, yloc_second_line],
             'k', clip_on=False)
    plt.plot([0, 1], [yloc_third_line, yloc_third_line],
             'k', clip_on=False)
    plt.plot([0, 1], [yloc_fourth_line, yloc_fourth_line],
             'k', clip_on=False)

    # Explanation of symbols: triangle is instrument, circle is mmi,
    # epicenter is black star
    # thick black line is rupture (if available)
    item_sep = [0.2, 0.28, 0.15]
    left_offset = 0.005
    label_pad = 0.02

    # Instrument
    triangle_marker_x = left_offset
    triangle_text_x = triangle_marker_x + label_pad
    plt.plot(triangle_marker_x, yloc_seventh_row, '^', markerfacecolor='w',
             markeredgecolor='k', markersize=6, mew=0.5, clip_on=False)
    plt.text(triangle_text_x,
             yloc_seventh_row,
             tdict['legend']['instrument'],
             va='center',
             ha='left')

    # Macroseismic
    circle_marker_x = triangle_text_x + item_sep[0]
    circle_text_x = circle_marker_x + label_pad
    plt.plot(circle_marker_x,
             yloc_seventh_row, 'o',
             markerfacecolor='w',
             markeredgecolor='k',
             markersize=4,
             mew=0.5)
    plt.text(circle_text_x,
             yloc_seventh_row,
             tdict['legend']['intensity'],
             va='center',
             ha='left')

    # Epicenter
    star_marker_x = circle_marker_x + item_sep[1]
    star_text_x = star_marker_x + label_pad
    plt.plot(star_marker_x,
             yloc_seventh_row, 'k*',
             markersize=12,
             mew=0.5)
    plt.text(star_text_x,
             yloc_seventh_row,
             tdict['legend']['epicenter'],
             va='center',
             ha='left')

    if not point_source:
        rup_marker_x = star_marker_x + item_sep[2]
        rup_text_x = rup_marker_x + label_pad
        rwidth = 0.02
        rheight = 0.05
        rup = patches.Rectangle(
            xy=(rup_marker_x - rwidth,
                yloc_seventh_row-0.5*rheight),
            width=rwidth,
            height=rheight,
            linewidth=2,
            edgecolor='k',
            facecolor='w'
        )
        cax.add_patch(rup)
        plt.text(rup_text_x,
                 yloc_seventh_row,
                 tdict['legend']['rupture'],
                 va='center',
                 ha='left')

    # Add conversion reference and shakemap version/process time
    version_x = 1.0
    tpl = (tdict['legend']['version'], map_version,
           tdict['legend']['processed'], process_time)
    plt.text(version_x, yloc_sixth_row,
             '%s %i: %s %s' % tpl,
             ha='right', va='center')

    ref = gmice.name
    refx = 0
    plt.text(refx, yloc_sixth_row,
             '%s %s' % (tdict['legend']['scale'], ref),
             va='center')

    nsteps = 10
    for i, width in enumerate(tdict['mmi_scale']['box_widths']):
        width /= 100
        textleft = sumwidth + width/2
        sumwidth += width
        plt.text(textleft, yloc_first_row,
                 tdict['mmi_scale']['shaking_labels'][i],
                 fontproperties=font1, **alignment)
        plt.text(textleft, yloc_second_row,
                 tdict['mmi_scale']['damage_labels'][i],
                 fontproperties=font1, **alignment)
        plt.text(textleft, yloc_third_row, acceleration[i],
                 fontproperties=font1, **alignment)
        plt.text(textleft, yloc_fourth_row, velocity[i],
                 fontproperties=font1, **alignment)

        if i == 0:
            font = font1
        else:
            font = font0
        th = plt.text(textleft, yloc_fifth_row,
                      tdict['mmi_scale']['intensity_labels'][i],
                      fontproperties=font, **alignment)
        th.set_path_effects([path_effects.Stroke(linewidth=2.0,
                                                 foreground='white'),
                             path_effects.Normal()])

        # draw right edge of cell
        plt.plot([gridleft+width, gridleft+width],
                 [bottom, top], 'k', clip_on=False)  # right

        # draw little colored rectangles inside the MMI cells
        if i > 0:
            left = gridleft
            ptop = yloc_fourth_line
            imt_min = imt_edges[i-1]
            imt_max = imt_edges[i]
            imts = np.linspace(imt_min, imt_max, nsteps)
            rights = np.linspace(gridleft, gridleft+width, nsteps)
            for mmi, right in zip(imts, rights):
                px = [left, right, right, left, left]
                py = [ptop, ptop, bottom, bottom, ptop]
                mmicolor = palette.getDataColor(mmi, color_format='hex')
                left = right
                plt.fill(px, py, mmicolor, ec=mmicolor)

        gridleft += width


def _draw_colorbar(fig, mmimap, tdict):
    """Draw an MMI colorbar in a separate axis from the map.

    Args:
        fig (Figure): Matplotlib Figure object.
        mmimap (ColorPalette): Impactutils MMI ColorPalette instance.
        tdict (dict): Dictionary containing the text strings in the user's
            choice of language.
    """
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

    plt.xticks(locs, tdict['mmi_scale']['mmi_colorbar_labels'])


def _label_close_to_edge(x, y, xmin, xmax, ymin, ymax):
    """Determine if a contour label is within specified distance of map edge.

    Args:
        x (float): Proposed X coordinate map label.
        y (float): Proposed Y coordinate map label.
        xmin (float): Left edge of map.
        xmax (float): Right edge of map.
        ymin (float): Bottom edge of map.
        ymax (float): Top edge of map.

    Returns:
        bool: True if label is "close" to edge of map.
    """
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


def _draw_graticules(ax, xmin, xmax, ymin, ymax):
    """Draw map graticules, tick labels on map axes.

    Args:
        ax (GeoAxes): Cartopy GeoAxes.
        xmin (float): Left edge of map (degrees).
        xmax (float): Right edge of map (degrees).
        ymin (float): Bottom edge of map (degrees).
        ymax (float): Bottom edge of map (degrees).
    """
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
    spans = {1: 0.25,
             2: 0.5,
             3: 1.0,
             5: 1.0,
             7: 2.0}

    span_keys = np.array(sorted(list(spans.keys())))

    nearest_xspan_idx = np.argmin(np.abs(int((xmax-xmin)) - span_keys))
    x_interval = spans[span_keys[nearest_xspan_idx]]

    nearest_yspan_idx = np.argmin(np.abs(int((ymax-ymin)) - span_keys))
    y_interval = spans[span_keys[nearest_yspan_idx]]

    # let's floor/ceil the edges to nearest 1/interval
    # gxmin = np.floor(xmin * interval) / interval
    # gxmax = np.ceil(xmax * interval) / interval
    # gymin = np.floor(ymin * interval) / interval
    # gymax = np.ceil(ymax * interval) / interval
    gxmin = x_interval*np.floor(xmin/x_interval)
    gxmax = x_interval*np.ceil(xmax/x_interval)
    gymin = y_interval*np.floor(ymin/y_interval)
    gymax = y_interval*np.ceil(ymax/y_interval)

    # check for meridian crossing
    crosses = False
    if gxmax < 0 and gxmax < gxmin:
        crosses = True
        gxmax += 360

    # shakemap way
    # ylocs = np.arange(np.floor(gymin), np.ceil(gymax) + interval, interval)
    # xlocs = np.arange(np.floor(gxmin), np.ceil(gxmax) + interval, interval)
    ylocs = np.arange(gymin, gymax+y_interval, y_interval)
    xlocs = np.arange(gxmin, gxmax+x_interval, x_interval)

    if crosses:
        xlocs[xlocs > 180] -= 360

    gl.xlocator = mticker.FixedLocator(xlocs)
    gl.ylocator = mticker.FixedLocator(ylocs)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}


def _get_shaded(ptopo, contour_colormap):
    """Get hill-shaded topo with topo colormap.

    Args:
        ptopo (ndarray): Projected topography numpy array.
        contour_colormap (ColorPalette): Topography color palette object.

    Returns:
        ndarray: Hill-shaded topography RGB image.
    """
    maxvalue = contour_colormap.vmax
    ls1 = LightSource(azdeg=300, altdeg=45)
    if np.allclose(ptopo, 0):
        intensity = np.full_like(ptopo, 1.0)
    else:
        maxvalue = contour_colormap.vmax
        ls2 = LightSource(azdeg=45, altdeg=45)
        intensity1 = ls1.hillshade(ptopo, fraction=0.25, vert_exag=VERT_EXAG)
        intensity2 = ls2.hillshade(ptopo, fraction=0.25, vert_exag=VERT_EXAG)
        intensity = intensity1 * 0.5 + intensity2 * 0.5

    ptoposc = ptopo / maxvalue
    rgba = contour_colormap.cmap(ptoposc)
    rgb = np.squeeze(rgba)

    draped_hsv = ls1.blend_hsv(rgb, np.expand_dims(intensity, 2))

    return draped_hsv


def _draw_title(imt, adict):
    """Draw the map title.
    Args:
        imt (str): IMT that is being drawn on the map ('MMI', 'PGV',
            'PGA', 'SA(x.y)').
        adict (dict): The dictionary containing the key geographic
            and ShakeMap data. See draw_map() for a descritption.
    """
    # Add a title
    tdict = adict['tdict']
    edict = adict['info']['input']['event_information']
    hlon = float(edict['longitude'])
    hlat = float(edict['latitude'])
    eloc = edict['location']
    try:
        etime = datetime.strptime(edict['origin_time'],
                                  constants.TIMEFMT)
    except ValueError:
        etime = datetime.strptime(edict['origin_time'],
                                  constants.ALT_TIMEFMT)
    timestr = etime.strftime(tdict['title_parts']['date_format'])
    mag = adict.get('display_magnitude')
    if mag is None:
        mag = float(edict['magnitude'])
    if hlon < 0:
        lonstr = '%s%.2f' % (tdict['title_parts']['west'], np.abs(hlon))
    else:
        lonstr = '%s%.2f' % (tdict['title_parts']['east'], hlon)
    if hlat < 0:
        latstr = '%s%.2f' % (tdict['title_parts']['south'], np.abs(hlat))
    else:
        latstr = '%s%.2f' % (tdict['title_parts']['north'], hlat)
    dep = float(edict['depth'])
    eid = edict['event_id']
    imtstr = tdict['IMTYPES'][imt]
    if len(eid) <= 10:
        fmt = ('%s\n%s %s: %s\n %s %s %s%.1f %s %s '
               '%s: %.1f%s %s:%s')
    else:
        fmt = ('%s\n%s %s: %s\n %s %s %s%.1f %s %s '
               '%s: %.1f%s\n%s:%s')
    tstr = fmt % (imtstr, adict['operator'],
                  tdict['title_parts']['shakemap'],
                  eloc, timestr,
                  tdict['title_parts']['timezone'],
                  tdict['title_parts']['magnitude'],
                  mag, latstr, lonstr,
                  tdict['title_parts']['depth'],
                  dep,
                  tdict['title_parts']['depth_units'],
                  tdict['title_parts']['event_id'],
                  eid)
    plt.title(tstr, fontsize=10, verticalalignment='bottom')


def _draw_stations(ax, stations, imt, intensity_colormap, geoproj, fill=True):
    """Draw station locations on the map.
    Args:
        ax (GeoAxes): Axes on which to draw stations.
        stations (dict): Dictionary containing station data.
        imt (str): IMT string.
        intensity_colormap (ColorPalette): MMI color palette object.
        geoproj (Proj): PlateCarree Pyproj Proj object.
        fill (bool): Whether or not to fill symbols.
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
            if np.isnan(mmi):
                mcolor = (mcolor[0], mcolor[1], mcolor[2], 0.0)
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
            if np.isnan(mmi):
                mcolor = (mcolor[0], mcolor[1], mcolor[2], 0.0)
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


def _get_draped(data, topodata, colormap):
    """Get array of data "draped" on topography.
    Args:
        data (ndarray): 2D Numpy array.
        topodata (ndarray): 2D Numpy array.
        colormap (ColorPalette): MMI color palette object.
    Returns:
        ndarray: Numpy array of data draped on topography.
    """

    maxvalue = colormap.vmax
    mmisc = data / maxvalue
    rgba_img = colormap.cmap(mmisc)
    rgb = np.squeeze(rgba_img[:, :, 0:3])

    if np.allclose(topodata, 0):
        intensity = np.full_like(topodata, 0.5)
    else:
        # use lightsource class to make our shaded topography
        ls1 = LightSource(azdeg=300, altdeg=45)
        ls2 = LightSource(azdeg=45, altdeg=45)
        intensity1 = ls1.hillshade(
            topodata, fraction=0.25, vert_exag=VERT_EXAG)
        intensity2 = ls2.hillshade(
            topodata, fraction=0.25, vert_exag=VERT_EXAG)
        intensity = intensity1 * 0.5 + intensity2 * 0.5

    ls = LightSource(azdeg=315, altdeg=45)
    draped_hsv = ls.blend_hsv(rgb, np.expand_dims(intensity, 2))

    return draped_hsv


def _clip_bounds(bbox, filename):
    """Clip input fiona-compatible vector file to input bounding box.

    Args:
        bbox (tuple): (xmin,ymin,xmax,ymax) desired clipping bounds.
        filename (str): Input name of file containing vector data in a
                        format compatible with fiona.
    Returns:
      BaseGeometry: Shapely Polygon or MultiPolygon.
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


def _draw_license(fig, adict):
    """Draw license information at the bottom of the figure if required.
    Args:
        fig (Figure): Matplotlib Figure object.
        adict (dict): The dictionary containing the key geographic
            and ShakeMap data. See draw_map() for a description.
    """
    logo_text = adict.get('license_text')
    if logo_text:
        lax = fig.add_axes([0.1, -0.05, 0.89, 0.04])
        logo_path = adict.get('license_logo')
        xpos = 0
        if logo_path:
            logo = image.imread(logo_path)
            h, w, colors = logo.shape
            ratio = w/h
            lax.imshow(logo, aspect='equal', extent=(0, ratio, 0, 1),
                       interpolation='bilinear')
            xpos = ratio + 0.25
        lax.set_aspect('equal', adjustable='box')
        lax.set_xlim(0, 20)
        lax.set_ylim(0, 1.025)
        lax.axis('off')
        from datetime import datetime
        year = datetime.now().strftime('%Y')
        text = logo_text.replace('%%YEAR%%', year)
        lax.text(xpos, 0.5, text, fontsize=9, va='center')


def draw_map(adict, override_scenario=False):
    """If adict['imtype'] is MMI, draw a map of intensity draped over
    topography, otherwise Draw IMT contour lines over hill-shaded topography.

    Args:
        adict (dictionary): A dictionary containing the following keys:
            'imtype' (str): The intensity measure type
            'topogrid' (Grid2d): A topography grid
            'allcities' (Cities): A list of global cities,
            'states_provinces' (Cartopy Feature): States/province boundaries.
            'countries' (Cartopy Feature): Country boundaries.
            'oceans' (Cartopy Feature): Oceans.
            'lakes' (Cartopy Feature): Lakes.
            'roads' (Shapely Feature): Roads.
            'faults' (Shapely Feature): Fault traces
            'datadir' (str): The path into which to deposit products
            'operator' (str): The producer of this shakemap
            'filter_size' (int): The size of the filter used before contouring
            'info' (dictionary): The shakemap info structure
            'component' (str): The intensity measure component being plotted
            'imtdict' (dictionary): Dict containing the IMT grids
            'rupdict' (dictionary): Dict containing the rupture data
            'stationdict' (dictionary): Dict of station data
            'config' (dictionary): The configuration data for this shakemap
            'tdict' (dictionary): The text strings to be printed on the map
                in the user's choice of language.
            'license_text' (str): License text to display at bottom of map
            'license_logo' (str): Path to license logo image to display
                next to license text
        override_scenario (bool): Turn off scenario watermark.

    Returns:
        Tuple of (Matplotlib figure, Matplotlib figure): Objects containing
        the map generated by this function, and the intensity legend,
        respectively. If the imtype of this map is not 'MMI', the second
        element of the tuple will be None.
    """
    imtype = adict['imtype']
    imtdict = adict['imtdict']      # mmidict
    imtdata = np.nan_to_num(imtdict['mean'], nan=0.0) # mmidata
    gd = GeoDict(imtdict['mean_metadata'])
    imtgrid = Grid2D(imtdata, gd)   # mmigrid

    gd = imtgrid.getGeoDict()

    # Retrieve the epicenter - this will get used on the map
    rupture = rupture_from_dict(adict['ruptdict'])
    origin = rupture.getOrigin()
    center_lat = origin.lat
    center_lon = origin.lon

    # load the cities data, limit to cities within shakemap bounds
    cities = adict['allcities'].limitByBounds((gd.xmin, gd.xmax,
                                               gd.ymin, gd.ymax))

    # get the map boundaries and figure size
    bounds, figsize, aspect = _get_map_info(gd)

    # Note: dimensions are: [left, bottom, width, height]
    dim_left = 0.1
    dim_bottom = 0.19
    dim_width = 0.8
    dim_height = dim_width/aspect
    if dim_height > 0.8:
        dim_height = 0.8
        dim_width = 0.8 * aspect
        dim_left = (1.0 - dim_width) / 2

    # Create the MercatorMap object, which holds a separate but identical
    # axes object used to determine collisions between city labels.
    mmap = MercatorMap(
        bounds, figsize, cities, padding=0.5,
        dimensions=[dim_left, dim_bottom, dim_width, dim_height])
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

    # get the projected IMT and topo grids
    pimtgrid, ptopogrid = _get_projected_grids(imtgrid, adict['topogrid'],
                                               projstr)

    # get the projected geodict
    proj_gd = pimtgrid.getGeoDict()

    pimtdata = pimtgrid.getData()
    ptopo_data = ptopogrid.getData()

    mmimap = ColorPalette.fromPreset('mmi')

    if imtype == 'MMI':
        draped_hsv = _get_draped(pimtdata, ptopo_data, mmimap)
    else:
        # get the draped topo data
        topo_colormap = ColorPalette.fromPreset('shaketopo')
        draped_hsv = _get_shaded(ptopo_data, topo_colormap)
        # convert units
        if imtype == 'PGV':
            pimtdata = np.exp(pimtdata)
        else:
            pimtdata = np.exp(pimtdata) * 100

    plt.sca(ax)
    ax.set_xlim(proj_gd.xmin, proj_gd.xmax)
    ax.set_ylim(proj_gd.ymin, proj_gd.ymax)
    img_extent = (proj_gd.xmin, proj_gd.xmax, proj_gd.ymin, proj_gd.ymax)

    plt.imshow(draped_hsv, origin='upper', extent=img_extent,
               zorder=IMG_ZORDER, interpolation='none')

    config = adict['config']
    gmice = get_object_from_config('gmice', 'modeling', config)
    gmice_imts = gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES
    gmice_pers = gmice.DEFINED_FOR_SA_PERIODS

    oqimt = imt.from_string(imtype)

    if imtype != 'MMI' and (not isinstance(oqimt, tuple(gmice_imts)) or
                            (isinstance(oqimt, imt.SA) and
                             oqimt.period not in gmice_pers)):
        my_gmice = None
    else:
        my_gmice = gmice

    if imtype != 'MMI':
        # call the contour module in plotting to get the vertices of the
        # contour lines
        contour_objects = contour(imtdict, imtype, adict['filter_size'],
                                  my_gmice)

        # get a color palette for the levels we have
        # levels = [c['properties']['value'] for c in contour_objects]

        # cartopy shapely feature has some weird behaviors, so I had to go
        # rogue and draw contour lines/labels myself.

        # To choose which contours to label, we will keep track of the lengths
        # of contours, grouped by isovalue
        contour_lens = defaultdict(lambda: [])
        def arclen(path):
            """
            Compute the arclength of *path*, which should be a list of pairs
            of numbers.
            """
            x0, y0 = [np.array(c) for c in zip(*path)]
            x1, y1 = [np.roll(c, -1) for c in (x0, y0)] # offset by 1
            # don't include first-last vertices as an edge:
            x0, y0, x1, y1 = [c[:-1] for c in (x0, y0, x1, y1)]
            return np.sum(np.sqrt((x0 - x1)**2 + (y0 - y1)**2))

        # draw dashed contours first, the ones over land will be overridden by
        # solid contours
        for contour_object in contour_objects:
            props = contour_object['properties']
            multi_lines = sShape(contour_object['geometry'])
            pmulti_lines = proj.project_geometry(multi_lines, src_crs=geoproj)
            for multi_line in pmulti_lines:
                pmulti_line = mapping(multi_line)['coordinates']
                x, y = zip(*pmulti_line)
                contour_lens[props['value']].append(arclen(pmulti_line))
                # color = imt_cmap.getDataColor(props['value'])
                ax.plot(x, y, color=props['color'], linestyle='dashed',
                        zorder=DASHED_CONTOUR_ZORDER)

        white_box = dict(
            boxstyle="round",
            ec=(0, 0, 0),
            fc=(1., 1, 1),
            color='k'
        )

        # draw solid contours next - the ones over water will be covered by
        # ocean polygon
        for contour_object in contour_objects:
            props = contour_object['properties']
            multi_lines = sShape(contour_object['geometry'])
            pmulti_lines = proj.project_geometry(multi_lines, src_crs=geoproj)

            # only label long contours (relative to others with the same
            # isovalue)
            min_len = np.array(contour_lens[props['value']]).mean()

            for multi_line in pmulti_lines:
                pmulti_line = mapping(multi_line)['coordinates']
                x, y = zip(*pmulti_line)
                # color = imt_cmap.getDataColor(props['value'])
                ax.plot(x, y, color=props['color'], linestyle='solid',
                        zorder=CONTOUR_ZORDER)
                if arclen(pmulti_line) >= min_len:
                    # try to label each segment with black text in a white box
                    xc = x[int(len(x)/3)]
                    yc = y[int(len(y)/3)]
                    if _label_close_to_edge(
                            xc, yc, proj_gd.xmin, proj_gd.xmax,
                            proj_gd.ymin, proj_gd.ymax):
                        continue
                    # TODO: figure out if box is going to go outside the map,
                    # if so choose a different point on the line.

                    # For small values, use scientific notation with 1 sig fig
                    # to avoid multiple contours labelled 0.0:
                    value = props['value']
                    fmt = '%.1g' if abs(value) < 0.1 else '%.1f'
                    ax.text(xc, yc, fmt % value, size=8,
                            ha="center", va="center",
                            bbox=white_box, zorder=AXES_ZORDER-1)

    # make the border thicker
    lw = 2.0
    ax.outline_patch.set_zorder(BORDER_ZORDER)
    ax.outline_patch.set_linewidth(lw)
    ax.outline_patch.set_joinstyle('round')
    ax.outline_patch.set_capstyle('round')

    # Coastlines will get drawn when we draw the ocean edges
    # ax.coastlines(resolution="10m", zorder=COAST_ZORDER, linewidth=3)

    if adict['states_provinces']:
        ax.add_feature(adict['states_provinces'], edgecolor='0.5',
                       zorder=COAST_ZORDER)

    if adict['countries']:
        ax.add_feature(adict['countries'], edgecolor='black',
                       zorder=BORDER_ZORDER)

    if adict['oceans']:
        ax.add_feature(adict['oceans'], edgecolor='black',
                       zorder=OCEAN_ZORDER)

    if adict['lakes']:
        ax.add_feature(adict['lakes'], edgecolor='black',
                       zorder=OCEAN_ZORDER)

    if adict['faults'] is not None:
        ax.add_feature(adict['faults'], edgecolor='firebrick',
                       zorder=ROAD_ZORDER)

    if adict['roads'] is not None:
        ax.add_feature(adict['roads'], edgecolor='dimgray',
                       zorder=ROAD_ZORDER)

    # draw graticules, ticks, tick labels
    _draw_graticules(ax, *bounds)

    # is this event a scenario?
    info = adict['info']
    etype = info['input']['event_information']['event_type']
    is_scenario = etype == 'SCENARIO'

    if is_scenario and not override_scenario:
        plt.text(
            center_lon, center_lat,
            adict['tdict']['title_parts']['scenario'],
            fontsize=72,
            zorder=SCENARIO_ZORDER, transform=geoproj,
            alpha=WATERMARK_ALPHA, color=WATERMARK_COLOR,
            horizontalalignment='center',
            verticalalignment='center',
            rotation=45,
            path_effects=[
                path_effects.Stroke(linewidth=1, foreground='black')]
        )

    # Draw the map scale in the unoccupied lower corner.
    corner = 'll'
    draw_scale(ax, corner, pady=0.05, padx=0.05, zorder=SCALE_ZORDER)

    # draw cities
    mmap.drawCities(shadow=True, zorder=CITIES_ZORDER, draw_dots=True)

    # Draw the epicenter as a black star
    plt.sca(ax)
    plt.plot(center_lon, center_lat, 'k*', markersize=16,
             zorder=EPICENTER_ZORDER, transform=geoproj)

    # draw the rupture polygon(s) in black, if not point rupture
    point_source = True
    if not isinstance(rupture, PointRupture):
        point_source = False
        json_dict = rupture._geojson
        for feature in json_dict['features']:
            for coords in feature['geometry']['coordinates']:
                for pcoords in coords:
                    poly2d = sLineString([xy[0:2] for xy in pcoords])
                    ppoly = proj.project_geometry(poly2d)
                    mppoly = mapping(ppoly)['coordinates']
                    for spoly in mppoly:
                        x, y = zip(*spoly)
                        ax.plot(x, y, 'k', lw=1, zorder=FAULT_ZORDER)

    # draw the station data on the map
    stations = adict['stationdict']
    _draw_stations(ax, stations, imtype, mmimap, geoproj)

    _draw_title(imtype, adict)

    process_time = info['processing']['shakemap_versions']['process_time']
    map_version = int(info['processing']['shakemap_versions']['map_version'])
    if imtype == 'MMI':
        _draw_mmi_legend(fig, mmimap, gmice, process_time,
                         map_version, point_source, adict['tdict'])
        # make a separate MMI legend
        fig2 = plt.figure(figsize=figsize)
        _draw_mmi_legend(fig2, mmimap, gmice, process_time,
                         map_version, point_source, adict['tdict'])

    else:
        _draw_imt_legend(fig, mmimap, imtype, gmice, process_time, map_version,
                         point_source, adict['tdict'])
        plt.draw()
        fig2 = None

    _draw_license(fig, adict)

    return (fig, fig2)
