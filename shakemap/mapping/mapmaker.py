# stdlib imports
from datetime import datetime
import os.path
from configobj import ConfigObj

# third party imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LightSource
import matplotlib.patheffects as path_effects
from matplotlib.font_manager import FontProperties
from matplotlib import patches

import cartopy.crs as ccrs  # projections
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import ShapelyFeature
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import pyproj

from shapely.geometry import shape as sShape
from shapely.geometry import Polygon as sPolygon
from shapely.geometry import GeometryCollection
from shapely.geometry import mapping

import fiona
from openquake.hazardlib import imt

# neic imports
from impactutils.mapping.mercatormap import MercatorMap
from impactutils.mapping.city import Cities
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
from shakelib.utils.imt_string import oq_to_file
from shakelib.gmice.wgrw12 import WGRW12
from shakemap.utils.utils import get_object_from_config
from shakemap.utils.config import get_config_paths

# define some constants
WATERCOLOR = '#7AA1DA'
FIGWIDTH = 10.0
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

IMTYPES = {
    'MMI': 'Macroseismic Intensity Map',
    'PGV': 'Peak Ground Velocity Map',
    'PGA': 'Peak Ground Acceleration Map',
    'SA(0.3)': '0.3 Second Peak Spectral Acceleration Map',
    'SA(1.0)': '1.0 Second Peak Spectral Acceleration Map',
    'SA(3.0)': '3.0 Second Peak Spectral Acceleration Map'
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

    simtgrid = imtgrid.interpolateToGrid(sampledict)

    # get topo layer and project it
    topogrid = topobase.interpolateToGrid(sampledict)

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

    # Map aspect
    aspect = pwidth/pheight

    fig_aspect = 1.0/(0.19 + 0.8/aspect)
    figheight = FIGWIDTH/fig_aspect
    bounds = (xmin, xmax, ymin, ymax)
    figsize = (FIGWIDTH, figheight)
    return (bounds, figsize, aspect)


def _draw_imt_legend(fig, palette, imtype, gmice):
    """Create a legend axis for non MMI plots.

    Args:
        fig (Figure): Matplotlib Figure object.
        levels (sequence): Sequence of contour levels.
        palette (ColorPalette): ColorPalette using range of input data and
            IMT_CMAP.
        imtype (str): One of 'PGV','PGA','SA(0.3)',etc.
    """
    units = {
        'PGV': '(cm/s)',
        'PGA': '(%g)',
        'SA(0.3)': '(%g)',
        'SA(1.0)': '(%g)',
        'SA(3.0)': '(%g)'
    }
    imtlabel = imtype + ' ' + units[imtype]
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


def _draw_mmi_legend(fig, palette, gmice, process_time, map_version,
                     point_source):
    """Create a legend axis for MMI plots.

    Args:
        fig (Figure): Matplotlib Figure object.
        palette (ColorPalette): ColorPalette using range of input data and
            IMT_CMAP.
        gmice: A gmice object.
        process_time (str): Process time.
        map_version (int): ShakeMap version.
        point_source (bool): Is the rupture a PointRupture?

    """
    cax = fig.add_axes([0.1, 0.00, 0.8, 0.15])
    plt.axis('off')
    cax_xmin, cax_xmax = cax.get_xlim()
    bottom, top = cax.get_ylim()
    plt.xlim(cax_xmin, cax_xmax)
    plt.ylim(bottom, top)

    shaking = [
        'SHAKING',
        'Not felt',
        'Weak',
        'Light',
        'Moderate',
        'Strong',
        'Very strong',
        'Severe',
        'Violent',
        'Extreme']
    damage = [
        'DAMAGE',
        'None',
        'None',
        'None',
        'Very light',
        'Light',
        'Moderate',
        'Moderate/heavy',
        'Heavy',
        'Very heavy']

    acceleration = ['PGA(%g)']
    velocity = ['PGV(cm/s)']

    imt_edges = np.array([0.5, 1.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
    intensities = [
        'INTENSITY', 'I', 'II-III',
        'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X+'
    ]

    widths = np.array([11.5, 7.75, 6.75, 7.0, 10.25,
                       8.5, 12.0, 16.25, 8.25, 11.75])/100

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
             'Seismic Instrument',
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
             'Macroseismic Observation',
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
             'Epicenter',
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
                 'Rupture',
                 va='center',
                 ha='left')

    # Add conversion reference and shakemap version/process time
    version_x = 1.0
    tpl = (map_version, process_time)
    plt.text(version_x, yloc_sixth_row,
             'Version %i: Processed %s' % tpl,
             ha='right', va='center')

    ref = gmice.name
    refx = 0
    plt.text(refx, yloc_sixth_row,
             'Scale based on %s' % ref,
             va='center')

    nsteps = 10
    for i in range(0, len(widths)):
        width = widths[i]
        textleft = sumwidth + width/2
        sumwidth += width
        plt.text(textleft, yloc_first_row, shaking[i],
                 fontproperties=font1, **alignment)
        plt.text(textleft, yloc_second_row, damage[i],
                 fontproperties=font1, **alignment)
        plt.text(textleft, yloc_third_row, acceleration[i],
                 fontproperties=font1, **alignment)
        plt.text(textleft, yloc_fourth_row, velocity[i],
                 fontproperties=font1, **alignment)

        if i == 0:
            font = font1
        else:
            font = font0
        th = plt.text(textleft, yloc_fifth_row, intensities[i],
                      fontproperties=font, **alignment)
        th.set_path_effects([path_effects.Stroke(linewidth=2.0,
                                                 foreground='white'),
                             path_effects.Normal()])

        # draw right edge of cell
        plt.plot([gridleft+widths[i], gridleft+widths[i]],
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

        gridleft += widths[i]


def _draw_colorbar(fig, mmimap):
    """Draw an MMI colorbar in a separate axis from the map.

    Args:
        fig (Figure): Matplotlib Figure object.
        mmimap (ColorPalette): Impactutils MMI ColorPalette instance.
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
    labels = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X']

    plt.xticks(locs, labels)


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
    ls2 = LightSource(azdeg=45, altdeg=45)
    intensity1 = ls1.hillshade(ptopo, fraction=0.25, vert_exag=VERT_EXAG)
    intensity2 = ls2.hillshade(ptopo, fraction=0.25, vert_exag=VERT_EXAG)
    intensity = intensity1 * 0.5 + intensity2 * 0.5

    ptoposc = ptopo / maxvalue
    rgba = contour_colormap.cmap(ptoposc)
    rgb = np.squeeze(rgba)

    draped_hsv = ls1.blend_hsv(rgb, np.expand_dims(intensity, 2))

    return draped_hsv


def _draw_title(imt, container, operator):
    """Draw the map title.
    Args:
        imt (str): IMT that is being drawn on the map ('MMI', 'PGV',
            'PGA', 'SA(x.y)').
        container (ShakeMapOutputContainer): HDF container of ShakeMap output.
        operator (str): Configured ShakeMap operator (NEIC, CISN, etc.)
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
    imtstr = IMTYPES[imt]
    if len(eid) <= 10:
        fmt = ('%s\n%s ShakeMap: %s\n %s UTC M%.1f %s %s '
               'Depth: %.1fkm ID:%s')
    else:
        fmt = ('%s\n%s ShakeMap: %s\n %s UTC M%.1f %s %s '
               'Depth: %.1fkm\nID:%s')
    tstr = fmt % (imtstr, operator, eloc, timestr, mag, latstr,
                  lonstr, dep, eid)
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
    # use lightsource class to make our shaded topography
    ls = LightSource(azdeg=315, altdeg=45)

    ls1 = LightSource(azdeg=300, altdeg=45)
    ls2 = LightSource(azdeg=45, altdeg=45)
    intensity1 = ls1.hillshade(
        topodata, fraction=0.25, vert_exag=VERT_EXAG)

    intensity2 = ls2.hillshade(
        topodata, fraction=0.25, vert_exag=VERT_EXAG)

    intensity = intensity1 * 0.5 + intensity2 * 0.5

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


def draw_intensity(container, topobase, oceanfile, outpath, operator,
                   borderfile=None, override_scenario=False):
    """Create a contour map showing MMI contours over greyscale population.

    Args:
        container (ShakeMapOutputContainer): HDF container of ShakeMap output.
        topobase (Grid2D): Topography grid to trim and project.
        oceanfile (str): Path to file containing ocean vector data in a
                         format compatible with fiona.
        outpath (str): Directory where output intensity.pdf and
                       intensity.jpg files will be made.
        operator (str): Configured ShakeMap operator (NEIC, CISN, etc.)
        borderfile (str): Shapefile containing country/state borders.
        override_scenario (bool): Turn off scenario watermark.

    Returns:
        str: Path to intensity PDF file.
        str: Path to intensity JPG file.
    """
    # get the geodict for the ShakeMap
    comp = container.getComponents('MMI')[0]
    imtdict = container.getIMTGrids('MMI', comp)
    mmidata = imtdict['mean']
    mmidict = imtdict['mean_metadata']
    gd = GeoDict(mmidict)
    mmigrid = Grid2D(mmidata, gd)

    # Retrieve the epicenter - this will get used on the map
    rupture = rupture_from_dict(container.getRuptureDict())
    origin = rupture.getOrigin()
    center_lat = origin.lat
    center_lon = origin.lon

    # load the cities data, limit to cities within shakemap bounds
    allcities = Cities.fromDefault()
    cities = allcities.limitByBounds((gd.xmin, gd.xmax, gd.ymin, gd.ymax))

    # get the map boundaries and figure size
    bounds, figsize, aspect = _get_map_info(gd)

    # Note: dimensions are: [left, bottom, width, height]
    dim_left = 0.1
    dim_bottom = 0.19
    dim_width = 0.8
    dim_height = dim_width/aspect

    # Create the MercatorMap object, which holds a separate but identical
    # axes object used to determine collisions between city labels.
    mmap = MercatorMap(
        bounds, figsize, cities, padding=0.5,
        dimensions=[dim_left, dim_bottom, dim_width, dim_height])
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

    # this is a workaround to an occasional problem where some vector layers
    # are not rendered. See
    # https://github.com/SciTools/cartopy/issues/1155#issuecomment-432941088
    proj._threshold /= 6

    # get the projected MMI and topo grids
    pmmigrid, ptopogrid = _get_projected_grids(mmigrid, topobase, projstr)

    # get the projected geodict
    proj_gd = pmmigrid.getGeoDict()

    # Use our GMT-inspired palette class to create population and MMI colormaps
    mmimap = ColorPalette.fromPreset('mmi')

    # drape the intensity data over the topography
    pmmi_data = pmmigrid.getData()
    ptopo_data = ptopogrid.getData()
    draped_hsv = _get_draped(pmmi_data, ptopo_data, mmimap)

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

    proj._threshold /= 6

    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
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

    _draw_stations(ax, stations, 'MMI', mmimap, geoproj)

    # Draw the epicenter as a black star
    plt.sca(ax)
    plt.plot(center_lon, center_lat, 'k*', markersize=16,
             zorder=EPICENTER_ZORDER, transform=geoproj)

    # draw the map title
    _draw_title('MMI', container, operator)

    # is this event a scenario?
    info = container.getMetadata()
    etype = info['input']['event_information']['event_type']
    is_scenario = etype == 'SCENARIO'

    if is_scenario and not override_scenario:
        plt.text(
            center_lon, center_lat, 'SCENARIO', fontsize=72,
            zorder=SCENARIO_ZORDER, transform=geoproj,
            alpha=WATERMARK_ALPHA, color=WATERMARK_COLOR,
            horizontalalignment='center',
            verticalalignment='center',
            rotation=45,
            path_effects=[path_effects.Stroke(
                linewidth=1,
                foreground='black')]
        )

    # draw the rupture polygon(s) in black, if not point rupture
    point_source = True
    if not isinstance(rupture, PointRupture):
        point_source = False
        json_dict = rupture._geojson
        for feature in json_dict['features']:
            rup_shape = sShape(feature['geometry'])
            sfeature = cfeature.ShapelyFeature(rup_shape, geoproj)
            ax.add_feature(sfeature, zorder=FAULT_ZORDER,
                           lw=1, edgecolor='k', facecolor=(0, 0, 0, 0))

    # draw graticules, ticks, tick labels
    _draw_graticules(ax, *bounds)

    # draw a separate intensity colorbar in a separate axes
    # _draw_colorbar(fig, mmimap)
    config = container.getConfig()
    gmice = get_object_from_config('gmice', 'modeling', config)
    process_time = info['processing']['shakemap_versions']['process_time']
    map_version = int(info['processing']['shakemap_versions']['map_version'])
    _draw_mmi_legend(fig, mmimap, gmice, process_time,
                     map_version, point_source)

    # make the map border thicker
    plt.sca(ax)
    lw = 2.0
    ax.outline_patch.set_zorder(BORDER_ZORDER)
    ax.outline_patch.set_linewidth(lw)
    ax.outline_patch.set_joinstyle('round')
    ax.outline_patch.set_capstyle('round')

    # ------------------------------------------ #
    # ***** Temp stuff for drawing circles ***** #
    # ------------------------------------------ #
    # Note this expects a file named 'circles.conf' to be located in the
    # event's 'current' directory. It can have a structure as follows, where
    # the radius is given in km:
    #
    # [line1]
    #     radius = 60
    #     marker = --k
    #     width = 2.0
    # [line2]
    #     radius = 120
    #     marker = --k
    #     width = 1.0
    #
    install_path, data_path = get_config_paths()
    datadir = os.path.join(
        data_path, info['input']['event_information']['id'], 'current')
    circle_conf = os.path.join(datadir, 'circles.conf')
    if os.path.isfile(circle_conf):
        cir_conf = ConfigObj(circle_conf)
        for k, v in cir_conf.items():
            # convert radius from km to m
            radius = float(v['radius']) * 1000
            pproj = pyproj.Proj(proj.proj4_init)
            cx, cy = pproj(origin.lon, origin.lat)
            npts = 500
            rad = np.linspace(0, 2*np.pi, npts)
            cir_x = np.cos(rad) * radius + cx
            cir_y = np.sin(rad) * radius + cy
            ax.plot(cir_x, cir_y, v['marker'], linewidth=float(v['width']))
    # ---------------------------------------------- #
    # ***** End temp stuff for drawing circles ***** #
    # ---------------------------------------------- #

    # create pdf and png output file names
    pdf_file = os.path.join(outpath, 'intensity.pdf')
    jpg_file = os.path.join(outpath, 'intensity.jpg')

    # save to pdf/jpeg
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.savefig(jpg_file, bbox_inches='tight')

    # make a separate MMI legend
    fig2 = plt.figure(figsize=figsize)
    _draw_mmi_legend(fig2, mmimap, gmice, process_time,
                     map_version, point_source)
    legend_file = os.path.join(outpath, 'mmi_legend.png')
    plt.savefig(legend_file, bbox_inches='tight')

    return (pdf_file, jpg_file, legend_file)


def draw_contour(container, imtype, topobase, oceanfile, outpath,
                 operator, filter_size, borderfile=None,
                 override_scenario=False):
    """Draw IMT contour lines over hill-shaded topography.

    Args:
        container (ShakeMapOutputContainer): HDF container of ShakeMap output.
        imtype (str): Type of IMT to be rendered.
        topobase (Grid2D): Topography grid to trim and project.
        oceanfile (str): Path to file containing ocean vector data in a
                         format compatible with fiona.
        outpath (str): Directory where output intensity.pdf and
                       intensity.jpg files will be made.
        operator (str): Configured ShakeMap operator (NEIC, CISN, etc.)
        filter_size (str): Smoothing filter size for contouring algorithm.
        borderfile (str): Shapefile containing country/state borders.
        override_scenario (bool): Turn off scenario watermark.

    Returns:
        str: Path to intensity PDF file.
        str: Path to intensity JPG file.
    """
    comp = container.getComponents(imtype)[0]
    imtdict = container.getIMTGrids(imtype, comp)
    imtdata = imtdict['mean']
    gd = GeoDict(imtdict['mean_metadata'])
    imtgrid = Grid2D(imtdata, gd)

    gd = imtgrid.getGeoDict()

    # Retrieve the epicenter - this will get used on the map
    rupture = rupture_from_dict(container.getRuptureDict())
    origin = rupture.getOrigin()
    center_lat = origin.lat

    # load the cities data, limit to cities within shakemap bounds
    allcities = Cities.fromDefault()
    cities = allcities.limitByBounds((gd.xmin, gd.xmax, gd.ymin, gd.ymax))

    # get the map boundaries and figure size
    bounds, figsize, aspect = _get_map_info(gd)

    # Note: dimensions are: [left, bottom, width, height]
    dim_left = 0.1
    dim_bottom = 0.19
    dim_width = 0.8
    dim_height = dim_width/aspect

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

    # this is a workaround to an occasional problem where some vector layers
    # are not rendered. See
    # https://github.com/SciTools/cartopy/issues/1155#issuecomment-432941088
    proj._threshold /= 6

    # get the projected MMI and topo grids
    pimtgrid, ptopogrid = _get_projected_grids(imtgrid, topobase, projstr)

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
    hillshade = _get_shaded(ptopogrid.getData(), topo_colormap)

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

    # get a color palette for the levels we have
    # levels = [c['properties']['value'] for c in contour_objects]

    # cartopy shapely feature has some weird behaviors, so I had to go rogue
    # and draw contour lines/labels myself.
    # draw dashed contours first, the ones over land will be overridden by
    # solid contours
    npoints = []
    for contour_object in contour_objects:
        props = contour_object['properties']
        multi_lines = sShape(contour_object['geometry'])
        pmulti_lines = proj.project_geometry(multi_lines, src_crs=geoproj)
        for multi_line in pmulti_lines:
            pmulti_line = mapping(multi_line)['coordinates']
            x, y = zip(*pmulti_line)
            npoints.append(len(x))
            # color = imt_cmap.getDataColor(props['value'])
            ax.plot(x, y, color=props['color'], linestyle='dashed',
                    zorder=DASHED_CONTOUR_ZORDER)

    white_box = dict(
        boxstyle="round",
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
            # color = imt_cmap.getDataColor(props['value'])
            ax.plot(x, y, color=props['color'], linestyle='solid',
                    zorder=CONTOUR_ZORDER)
            if len(x) > min_npoints:
                # try to label each segment with black text in a white box
                xc = x[int(len(x)/3)]
                yc = y[int(len(y)/3)]
                if _label_close_to_edge(
                        xc, yc, proj_gd.xmin, proj_gd.xmax,
                        proj_gd.ymin, proj_gd.ymax):
                    continue
                # TODO: figure out if box is going to go outside the map, if so
                # choose a different point on the line.
                ax.text(xc, yc, '%.1f' % props['value'], size=8,
                        ha="center", va="center",
                        bbox=white_box, zorder=AXES_ZORDER-1)

    # make the border thicker
    lw = 2.0
    ax.outline_patch.set_zorder(BORDER_ZORDER)
    ax.outline_patch.set_linewidth(lw)
    ax.outline_patch.set_joinstyle('round')
    ax.outline_patch.set_capstyle('round')

    # clip the ocean data to the shakemap
    bbox = (gd.xmin, gd.ymin, gd.xmax, gd.ymax)
    oceanshapes = _clip_bounds(bbox, oceanfile)

    ax.add_feature(
        ShapelyFeature(oceanshapes, crs=geoproj),
        facecolor=WATERCOLOR,
        zorder=OCEAN_ZORDER
    )

    # draw 10m res coastlines
    ax.coastlines(resolution="10m", zorder=COAST_ZORDER, linewidth=3)

    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

    ax.add_feature(states_provinces, edgecolor='black', zorder=COAST_ZORDER)

    # draw graticules, ticks, tick labels
    _draw_graticules(ax, *bounds)

    # is this event a scenario?
    info = container.getMetadata()
    etype = info['input']['event_information']['event_type']
    is_scenario = etype == 'SCENARIO'

    # Retrieve the epicenter - this will get used on the map
    center_lat = origin.lat
    center_lon = origin.lon

    if is_scenario and not override_scenario:
        plt.text(
            center_lon, center_lat, 'SCENARIO', fontsize=72,
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
    mmap.drawCities(
        shadow=True, zorder=CITIES_ZORDER, draw_dots=True)

    # Draw the epicenter as a black star
    plt.sca(ax)
    plt.plot(center_lon, center_lat, 'k*', markersize=16,
             zorder=EPICENTER_ZORDER, transform=geoproj)

    # draw the rupture polygon(s) in black, if not point rupture
    rupture = rupture_from_dict(container.getRuptureDict())
    if not isinstance(rupture, PointRupture):
        json_dict = rupture._geojson
        # shapes = []
        for feature in json_dict['features']:
            rup_shape = sShape(feature['geometry'])
            sfeature = cfeature.ShapelyFeature(rup_shape, geoproj)
            ax.add_feature(sfeature, zorder=FAULT_ZORDER,
                           lw=1, edgecolor='k', facecolor=(0, 0, 0, 0))

    # draw the station data on the map
    stations = container.getStationDict()
    mmimap = ColorPalette.fromPreset('mmi')
    _draw_stations(ax, stations, imtype, mmimap, geoproj)

    _draw_title(imtype, container, operator)

    _draw_imt_legend(fig, mmimap, imtype, gmice)

    # ------------------------------------------ #
    # ***** Temp stuff for drawing circles ***** #
    # ------------------------------------------ #
    # Note this expects a file named 'circles.conf' to be located in the
    # event's 'current' directory. It can have a structure as follows, where
    # the radius is given in km:
    #
    # [line1]
    #     radius = 60
    #     marker = --k
    #     width = 2.0
    # [line2]
    #     radius = 120
    #     marker = --k
    #     width = 1.0
    #
    install_path, data_path = get_config_paths()
    datadir = os.path.join(
        data_path, info['input']['event_information']['id'], 'current')
    circle_conf = os.path.join(datadir, 'circles.conf')
    if os.path.isfile(circle_conf):
        cir_conf = ConfigObj(circle_conf)
        for k, v in cir_conf.items():
            # convert radius from km to m
            radius = float(v['radius']) * 1000
            pproj = pyproj.Proj(proj.proj4_init)
            cx, cy = pproj(origin.lon, origin.lat)
            npts = 500
            rad = np.linspace(0, 2*np.pi, npts)
            cir_x = np.cos(rad) * radius + cx
            cir_y = np.sin(rad) * radius + cy
            ax.plot(cir_x, cir_y, v['marker'], linewidth=float(v['width']))
    # ---------------------------------------------- #
    # ***** End temp stuff for drawing circles ***** #
    # ---------------------------------------------- #

    # save plot to file
    fileimt = oq_to_file(imtype)
    plt.draw()
    pdf_file = os.path.join(outpath, '%s.pdf' % (fileimt))
    jpg_file = os.path.join(outpath, '%s.jpg' % (fileimt))
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.savefig(jpg_file, bbox_inches='tight')

    return (pdf_file, jpg_file)
