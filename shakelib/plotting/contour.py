# usgs imports
from impactutils.colors.cpalette import ColorPalette

# third party imports
from shapely.geometry import MultiLineString
from shapely.geometry import mapping
from scipy.ndimage.filters import median_filter
from skimage import measure
import numpy as np


def contour(container, imtype, component, filter_size):
    """
    Generate contours of a specific IMT and return as a Shapely
    MultiLineString object.

    Args:
        container (ShakeMapOutputContainer): ShakeMapOutputContainer
            with ShakeMap output data.
        imtype (str): String containing the name of an Intensity
            Measure Type found in container.
        component (str): Intensity Measure component found in container.
        filter_size (int): Integer filter (see
            https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.ndimage.filters.median_filter.html)
    Returns:
        list: List of dictionaries containing two fields

                - geometry: GeoJSON-like representation of one of the objects
                  in https://toblerity.org/fiona/manual.html#geometry-types
                - properties: Dictionary of properties describing that
                  feature.

    Raises:
        NotImplementedError -- if the user attempts to contour a data file
            with sets of points rather than grids.
    """  # noqa
    intensity_colormap = ColorPalette.fromPreset('mmi')
    imtdict = container.getIMTGrids(imtype, component)
    gridobj = imtdict['mean']
    grid = gridobj.getData()
    metadata = gridobj.getGeoDict().asDict()
    if imtype == 'MMI':
        sgrid = grid
        units = 'mmi'
    elif imtype == 'PGV':
        sgrid = np.exp(grid)
        units = 'cms'
    else:
        sgrid = np.exp(grid) * 100.0
        units = 'pctg'
    if filter_size > 0:
        fgrid = median_filter(sgrid, size=filter_size)
    else:
        fgrid = sgrid

    interval_type = 'log'
    if imtype == 'MMI':
        interval_type = 'linear'
    intervals = getContourLevels(
        np.min(fgrid), np.max(fgrid), itype=interval_type)

    lonstart = metadata['xmin']
    latstart = metadata['ymin']
    lonspan = np.abs(metadata['xmax'] - lonstart)
    latspan = np.abs(metadata['ymax'] - latstart)
    nlon = metadata['nx']
    nlat = metadata['ny']

    line_strings = []  # dictionary of MultiLineStrings and props

    for cval in intervals:
        contours = measure.find_contours(fgrid, cval)
        #
        # Convert coords to geographic coordinates; the coordinates
        # are returned in row, column order (i.e., (y, x))
        #
        new_contours = []
        plot_contours = []
        for ic, coords in enumerate(contours):  # coords is a line segment
            if len(coords) <= 20:  # skipping little contour islands?
                continue

            mylons = coords[:, 1] * lonspan / nlon + lonstart
            mylats = (nlat - coords[:, 0]) * latspan / nlat + latstart
            contours[ic][:, 0] = mylons[:]
            contours[ic][:, 1] = mylats[:]
            plot_contours.append(contours[ic])
            new_contours.append(contours[ic].tolist())

        if len(new_contours):
            mls = MultiLineString(new_contours)
            props = {
                'value': cval,
                'units': units
            }
            if imtype == 'MMI':
                color_array = np.array(intensity_colormap.getDataColor(cval))
                color_rgb = np.array(
                    color_array[0:3] * 255, dtype=int).tolist()
                props['color'] = '#%02x%02x%02x' % tuple(color_rgb)
                if (cval * 2) % 2 == 1:
                    props['weight'] = 4
                else:
                    props['weight'] = 2
            line_strings.append(
                {
                    'geometry': mapping(mls),
                    'properties': props
                }
            )
    return line_strings


def getContourLevels(dmin, dmax, itype='log'):
    """
    Get contour levels given min/max values and desired IMT.

    Use itype='log' for any IMT that is logarithmically distributed, such as
    PGA, PGV, and Sa. Linear for MMI.

    Args:
        dmin (float): Minimum value of data to contour.
        dmax (float): Maximum value of data to contour.
        itype (str): Interval type; default is 'log', anythign else
            indicates linear intervals.

    Returns:
        ndarray: Numpy array of contour levels.

    """
    if itype == 'log':
        # Within-decade label values
        dec_inc = np.array([1, 2, 5], dtype=float)

        # Upper and lower decades
        lower_dec = np.floor(np.log10(dmin))
        upper_dec = np.ceil(np.log10(dmax))

        # Array of decades
        decades = np.arange(lower_dec, upper_dec + 1)

        # Construct levels
        levels = np.concatenate([np.power(10, d) * dec_inc for d in decades])
        levels = levels[(levels < dmax) & (levels > dmin)]
    else:
        # MMI contours are every 0.5 units
        levels = np.arange(
            np.ceil(dmin * 2) / 2,
            np.floor(dmax * 2) / 2 + 0.5,
            0.5
        )
    return levels
