import numpy as np

from openquake.hazardlib.geo.utils import OrthographicProjection

from shakelib.rupture.edge_rupture import EdgeRupture
from shakelib.rupture.quad_rupture import QuadRupture
from shakelib.rupture.base import Rupture

from strec.gmreg import Regionalizer

DEFAULT_ACTIVE_COEFFS = [27.24, 250.4, 579.1]
DEFAULT_STABLE_COEFFS = [63.4, 465.4, 581.3]


def get_extent(rupture, config=None):
    """
    Method to compute map extent from rupture.

    Args:
        rupture (Rupture): A ShakeMap Rupture instance.

    Returns:
        tuple: lonmin, lonmax, latmin, latmax rounded to the nearest
        arc-minute..

    """

    # check to see what parameters are specified in the extent config
    coeffs = []
    spans = {}
    bounds = []
    if config is not None:
        if 'extent' in config:
            if 'coefficients' in config['extent']:
                if 'coeffs' in config['extent']['coefficients']:
                    if config['extent']['coefficients']['coeffs'][0] != 0.0:
                        coeffs = config['extent']['coefficients']['coeffs']
            if 'magnitude_spans' in config['extent']:
                if len(config['extent']['magnitude_spans']):
                    if isinstance(config['extent']['magnitude_spans'], dict):
                        spans = config['extent']['magnitude_spans']
            if 'bounds' in config['extent']:
                if 'extent' in config['extent']['bounds']:
                    if config['extent']['bounds']['extent'][0] != -999.0:
                        bounds = config['extent']['bounds']['extent']

    if len(bounds):
        xmin, ymin, xmax, ymax = bounds
        return (xmin, xmax, ymin, ymax)

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

    if len(spans):
        xmin = None
        xmax = None
        ymin = None
        ymax = None
        for spankey, span in spans.items():
            if mag > span[0] and mag <= span[1]:
                ymin = clat - span[2]/2
                ymax = clat + span[2]/2
                xmin = clon - span[3]/2
                xmax = clon + span[3]/2
                break
        if xmin is not None:
            return (xmin, xmax, ymin, ymax)

    # Is this a stable or active tectonic event?
    # (this could be made an attribute of the ShakeMap Origin class)
    hypo = origin.getHypo()
    stable = is_stable(hypo.longitude, hypo.latitude)

    if stable is False:
        if mag < 6.48:
            mindist_km = 100.
        else:
            if len(coeffs):
                C1, C2, C3 = coeffs
            else:
                C1, C2, C3 = DEFAULT_ACTIVE_COEFFS
            mindist_km = C1 * mag**2 - C2 * mag + C3
    else:
        if mag < 6.10:
            mindist_km = 100.
        else:
            if len(coeffs):
                C1, C2, C3 = coeffs
            else:
                C1, C2, C3 = DEFAULT_STABLE_COEFFS
            mindist_km = C1 * mag**2 - C2 * mag + C3

    # Apply an upper limit on extent. This should only matter for large
    # magnitudes (> ~8.25) in stable tectonic environments.
    if mindist_km > 1000.:
        mindist_km = 1000.

    # Projection
    proj = OrthographicProjection(clon - 4, clon + 4, clat + 4, clat - 4)
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
    reg = Regionalizer.load()
    region_info = reg.getRegions(lat, lon, 0)
    return region_info['TectonicRegion'] == 'Stable'
