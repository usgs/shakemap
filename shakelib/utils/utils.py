import os
import json
import pkg_resources

import numpy as np
from shapely.geometry import Polygon, Point

from openquake.hazardlib.geo.utils import get_orthographic_projection

from shakelib.rupture.edge_rupture import EdgeRupture
from shakelib.rupture.quad_rupture import QuadRupture
from shakelib.rupture.base import Rupture


def get_extent(rupture):
    """
    Method to compute map extent from rupture.

    Args:
        rupture (Rupture): A ShakeMap Rupture instance.

    Returns:
        tuple: lonmin, lonmax, latmin, latmax rounded to the nearest
        arc-minute..

    """

    if not rupture or not isinstance(rupture, Rupture):
        raise TypeError('get_extent() takes exactly 1 argument (0 given)')

    origin = rupture.getOrigin()
    if isinstance(rupture, (QuadRupture, EdgeRupture)):
        lats = rupture.lats
        lons = rupture.lons

        # Remove nans
        lons = lons[~np.isnan(lons)]
        lats = lats[~np.isnan(lats)]

        clat = 0.5 * (np.nanmax(lats) + np.nanmin(lats))
        clon = 0.5 * (np.nanmax(lons) + np.nanmin(lons))
    else:
        clat = origin.lat
        clon = origin.lon

    mag = origin.mag

    # Is this a stable or active tectonic event?
    # (this could be made an attribute of the ShakeMap Origin class)
    hypo = origin.getHypo()
    stable = is_stable(hypo.longitude, hypo.latitude)

    if stable is False:
        if mag < 6.48:
            mindist_km = 100.
        else:
            mindist_km = 27.24 * mag**2 - 250.4 * mag + 579.1
    else:
        if mag < 6.10:
            mindist_km = 100.
        else:
            mindist_km = 63.4 * mag**2 - 465.4 * mag + 581.3

    # Apply an upper limit on extent. This should only matter for large
    # magnitudes (> ~8.25) in stable tectonic environments.
    if mindist_km > 1000.:
        mindist_km = 1000.

    # Projection
    proj = get_orthographic_projection(clon - 4, clon + 4, clat + 4, clat - 4)
    if isinstance(rupture, (QuadRupture, EdgeRupture)):
        ruptx, rupty = proj(lons, lats)
    else:
        ruptx, rupty = proj(clon, clat)

    xmin = np.nanmin(ruptx) - mindist_km
    ymin = np.nanmin(rupty) - mindist_km
    xmax = np.nanmax(ruptx) + mindist_km
    ymax = np.nanmax(rupty) + mindist_km

    # Put a limit on range of aspect ratio
    dx = xmax - xmin
    dy = ymax - ymin
    ar = dy / dx
    if ar > 1.25:
        # Inflate x
        dx_target = dy / 1.25
        ddx = dx_target - dx
        xmax = xmax + ddx / 2
        xmin = xmin - ddx / 2
    if ar < 0.6:
        # inflate y
        dy_target = dx * 0.6
        ddy = dy_target - dy
        ymax = ymax + ddy / 2
        ymin = ymin - ddy / 2

    lonmin, latmin = proj(np.array([xmin]), np.array([ymin]), reverse=True)
    lonmax, latmax = proj(np.array([xmax]), np.array([ymax]), reverse=True)

    #
    # Round coordinates to the nearest minute -- that should make the
    # output grid register with common grid resolutions (60c, 30c,
    # 15c, 7.5c)
    #
    return _round_coord(lonmin[0]), _round_coord(lonmax[0]), \
           _round_coord(latmin[0]), _round_coord(latmax[0])


def _round_coord(coord):
    """
    Round a number to the nearest arc-minute
    """
    dm = 1.0 / 60
    mm = coord / dm
    imm = int(mm + 0.5)
    return imm * dm


def is_stable(lon, lat):
    """
    Determine if point is located in the US stable tectonic region. Uses the
    same boundary as the US NSHMP and so this function needs to be modified to
    work outside of the US.

    Args:
        lon (float): Lognitude.
        lat (float): Latitude.

    Returns:
        bool: Is the point classified as tectonically stable.

    """
    p = Point((lon, lat))
    pfile = pkg_resources.resource_filename('shakelib.utils', 
            os.path.join('data', 'cratons.geojson'))
    with open(pfile) as f:
        cratons = json.load(f)
    coord_list = cratons['features'][0]['geometry']['coordinates']
    for clist in coord_list:
        poly = Polygon(clist[0])
        if p.within(poly):
            return True
    return False
