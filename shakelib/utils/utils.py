import numpy as np
import logging

from openquake.hazardlib import imt
from openquake.hazardlib import const
from openquake.hazardlib.geo.utils import OrthographicProjection
from openquake.hazardlib.gsim.base import SitesContext
from openquake.hazardlib.gsim.base import DistancesContext
from openquake.hazardlib.gsim.base import RuptureContext

from shakemap.utils.utils import get_object_from_config
from shakelib.rupture.edge_rupture import EdgeRupture
from shakelib.rupture.quad_rupture import QuadRupture
from shakelib.rupture.base import Rupture
from shakelib.multigmpe import MultiGMPE
from shakelib.sites import Sites
from shakelib.station import StationList

from strec.gmreg import Regionalizer

# ACR GMPE/GMICE
from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014
from openquake.hazardlib.gsim.campbell_bozorgnia_2014 import (
    CampbellBozorgnia2014)
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
from shakelib.gmice.wgrw12 import WGRW12

# SCR GMPE/GMICE
from openquake.hazardlib.gsim.frankel_1996 import FrankelEtAl1996MwNSHMP2008
from openquake.hazardlib.gsim.toro_1997 import ToroEtAl1997MwNSHMP2008
from openquake.hazardlib.gsim.silva_2002 import SilvaEtAl2002MwNSHMP2008
from openquake.hazardlib.gsim.campbell_2003 import Campbell2003MwNSHMP2008
from openquake.hazardlib.gsim.tavakoli_pezeshk_2005 import (
    TavakoliPezeshk2005MwNSHMP2008)
from openquake.hazardlib.gsim.atkinson_boore_2006 import (
    AtkinsonBoore2006Modified2011)
from openquake.hazardlib.gsim.pezeshk_2011 import PezeshkEtAl2011
from openquake.hazardlib.gsim.boore_atkinson_2011 import Atkinson2008prime
from openquake.hazardlib.gsim.somerville_2001 import (
    SomervilleEtAl2001NSHMP2008)
from shakelib.gmice.ak07 import AK07


DEFAULT_ACTIVE_COEFFS = [27.24, 250.4, 579.1]
DEFAULT_STABLE_COEFFS = [63.4, 465.4, 581.3]


def replace_dyfi(stationfile, dyfi_xml):
    """Remove any non-instrumented data from station file, add DYFI.

    Args:
        stationfile (str): Existing station data file, presumed to
                           contain old DYFI data.
        dyfi_xml (str): DYFI XML data file, which will be added to
                        instrumented station data.
    Returns:
        StationList: Object containing merged data.

    """
    stations = StationList.loadFromFiles([stationfile])
    # reach into the internal database and find the instrumented stations
    conn = stations.db
    cursor = stations.cursor
    query1 = ('SELECT id from station WHERE instrumented = '
              '0 and (network = "DYFI" or network="CIIM")')
    cursor.execute(query1)
    rows = cursor.fetchall()
    for row in rows:
        sid = row[0]
        query2 = 'DELETE FROM amp WHERE station_id="%s"' % sid
        cursor.execute(query2)
        conn.commit()
        query3 = 'DELETE FROM station where id="%s"' % sid
        cursor.execute(query3)
        conn.commit()

    # now insert the dyfi data
    stations.addData([dyfi_xml])
    return stations


def get_extent(rupture=None, config=None):
    """
    Method to compute map extent from rupture. There are numerous methods for
    getting the extent:
        - It can be specified directly in the config file,
        - it can be hard coded for specific magnitude ranges in the config
          file, or
        - it can be based on the MultiGMPE for the event.

    All methods except for the first requires a rupture object.

    If no config is provided then a rupture is required and the extent is based
    on a generic set of active/stable.

    Args:
        rupture (Rupture): A ShakeMap Rupture instance.
        config (ConfigObj): ShakeMap config object.

    Returns:
        tuple: lonmin, lonmax, latmin, latmax rounded to the nearest
        arc-minute..

    """

    # -------------------------------------------------------------------------
    # Check to see what parameters are specified in the extent config
    # -------------------------------------------------------------------------
    spans = {}
    bounds = []
    offsets = None
    if config is not None:
        if 'extent' in config:
            if 'magnitude_spans' in config['extent']:
                if len(config['extent']['magnitude_spans']):
                    if isinstance(config['extent']['magnitude_spans'], dict):
                        spans = config['extent']['magnitude_spans']
            if 'bounds' in config['extent']:
                if 'extent' in config['extent']['bounds']:
                    if config['extent']['bounds']['extent'][0] != -999.0:
                        bounds = config['extent']['bounds']['extent']
            if 'relative_offset' in config['extent']:
                if isinstance(config['extent']['relative_offset'], list):
                    offsets = config['extent']['relative_offset']

    # -------------------------------------------------------------------------
    # Simplest option: extent was specified in the config, use that and exit.
    # -------------------------------------------------------------------------
    if len(bounds):
        xmin, ymin, xmax, ymax = bounds
        return (xmin, xmax, ymin, ymax)

    if not rupture or not isinstance(rupture, Rupture):
        raise TypeError('get_extent() requires a rupture object if the extent '
                        'is not specified in the config object.')

    # -------------------------------------------------------------------------
    # Second simplest option: spans are hardcoded based on magnitude
    # -------------------------------------------------------------------------
    extent = None
    if len(spans):
        extent = _get_extent_from_spans(rupture, spans)

    # -------------------------------------------------------------------------
    # Otherwise, use MultiGMPE to get spans
    # -------------------------------------------------------------------------
    if extent is None:
        extent = _get_extent_from_multigmpe(rupture, config)

    if offsets is None:
        return extent

    # -------------------------------------------------------------------------
    # Apply relative offsets
    # -------------------------------------------------------------------------
    (xmin, xmax, ymin, ymax) = extent

    xspan = xmax - xmin
    yspan = ymax - ymin

    xmin += xspan * offsets[0]
    xmax += xspan * offsets[0]
    ymin += yspan * offsets[1]
    ymax += yspan * offsets[1]

    return (xmin, xmax, ymin, ymax)


def _get_extent_from_spans(rupture, spans=[]):
    """
    Choose extent based on magnitude using a hardcoded list of spans
    based on magnitude ranges.
    """
    (clon, clat) = _rupture_center(rupture)
    mag = rupture.getOrigin().mag
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
    return None


def _get_extent_from_multigmpe(rupture, config=None):
    """
    Use MultiGMPE to determine extent
    """
    (clon, clat) = _rupture_center(rupture)
    origin = rupture.getOrigin()
    if config is not None:
        gmpe = MultiGMPE.from_config(config)
        gmice = get_object_from_config('gmice', 'modeling', config)
        if imt.SA in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES:
            default_imt = imt.SA(1.0)
        elif imt.PGV in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES:
            default_imt = imt.PGV()
        else:
            default_imt = imt.PGA()
    else:
        # Put in some default values for conf
        config = {
            'extent': {
                'mmi': {
                    'threshold': 4.5,
                    'mindist': 100,
                    'maxdist': 1000
                }
            }
        }

        # Generic GMPEs choices based only on active vs stable
        # as defaults...
        stable = is_stable(origin.lon, origin.lat)
        if not stable:
            ASK14 = AbrahamsonEtAl2014()
            CB14 = CampbellBozorgnia2014()
            CY14 = ChiouYoungs2014()
            gmpes = [ASK14, CB14, CY14]
            site_gmpes = None
            weights = [1/3.0, 1/3.0, 1/3.0]
            gmice = WGRW12()
        else:
            Fea96 = FrankelEtAl1996MwNSHMP2008()
            Tea97 = ToroEtAl1997MwNSHMP2008()
            Sea02 = SilvaEtAl2002MwNSHMP2008()
            C03 = Campbell2003MwNSHMP2008()
            TP05 = TavakoliPezeshk2005MwNSHMP2008()
            AB06p = AtkinsonBoore2006Modified2011()
            Pea11 = PezeshkEtAl2011()
            Atk08p = Atkinson2008prime()
            Sea01 = SomervilleEtAl2001NSHMP2008()
            gmpes = [Fea96, Tea97, Sea02, C03,
                     TP05, AB06p, Pea11, Atk08p, Sea01]
            site_gmpes = [AB06p]
            weights = [0.16, 0.0, 0.0, 0.17, 0.17, 0.3, 0.2, 0.0, 0.0]
            gmice = AK07()

        gmpe = MultiGMPE.from_list(
            gmpes, weights, default_gmpes_for_site=site_gmpes)
        default_imt = imt.SA(1.0)

    min_mmi = config['extent']['mmi']['threshold']
    sd_types = [const.StdDev.TOTAL]

    # Distance context
    dx = DistancesContext()
    # This imposes minimum/ maximum distances of:
    #   80 and 800 km; could make this configurable
    d_min = config['extent']['mmi']['mindist']
    d_max = config['extent']['mmi']['maxdist']
    dx.rjb = np.logspace(np.log10(d_min), np.log10(d_max), 2000)
    # Details don't matter for this; assuming vertical surface rupturing fault
    # with epicenter at the surface.
    dx.rrup = dx.rjb
    dx.rhypo = dx.rjb
    dx.repi = dx.rjb
    dx.rx = np.zeros_like(dx.rjb)
    dx.ry0 = np.zeros_like(dx.rjb)
    dx.rvolc = np.zeros_like(dx.rjb)

    # Sites context
    sx = SitesContext()
    # Set to soft soil conditions
    sx.vs30 = np.full_like(dx.rjb, 180)
    sx = MultiGMPE.set_sites_depth_parameters(sx, gmpe)
    sx.vs30measured = np.full_like(sx.vs30, False, dtype=bool)
    sx = Sites._addDepthParameters(sx)
    sx.backarc = np.full_like(sx.vs30, False, dtype=bool)

    # Rupture context
    rx = RuptureContext()
    rx.mag = origin.mag
    rx.rake = 0.0
    # From WC94...
    rx.width = 10**(-0.76 + 0.27*rx.mag)
    rx.dip = 90.0
    rx.ztor = origin.depth
    rx.hypo_depth = origin.depth

    gmpe_imt_mean, _ = gmpe.get_mean_and_stddevs(
        sx, rx, dx, default_imt, sd_types)

    # Convert to MMI
    gmpe_to_mmi, _ = gmice.getMIfromGM(gmpe_imt_mean, default_imt)

    # Minimum distance that exceeds threshold MMI?
    dists_exceed_mmi = dx.rjb[gmpe_to_mmi > min_mmi]
    if len(dists_exceed_mmi):
        mindist_km = np.max(dists_exceed_mmi)
    else:
        mindist_km = d_min

    # Get a projection
    proj = OrthographicProjection(clon - 4, clon + 4, clat + 4, clat - 4)
    if isinstance(rupture, (QuadRupture, EdgeRupture)):
        ruptx, rupty = proj(
            rupture.lons[~np.isnan(rupture.lons)],
            rupture.lats[~np.isnan(rupture.lats)]
        )
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
    if ar > 1.2:
        # Inflate x
        dx_target = dy / 1.2
        ddx = dx_target - dx
        xmax = xmax + ddx / 2
        xmin = xmin - ddx / 2
    if ar < 0.83:
        # Inflate y
        dy_target = dx * 0.83
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
    logging.debug("Extent: %f, %f, %f, %f" %
                  (lonmin, lonmax, latmin, latmax))
    return _round_coord(lonmin[0]), _round_coord(lonmax[0]), \
        _round_coord(latmin[0]), _round_coord(latmax[0])


def _rupture_center(rupture):
    """
    Find the central point of a rupture
    """
    origin = rupture.getOrigin()
    if isinstance(rupture, (QuadRupture, EdgeRupture)):
        # For an extended rupture, it is the midpoint between the extent of the
        # verticies
        lats = rupture.lats
        lons = rupture.lons

        # Remove nans
        lons = lons[~np.isnan(lons)]
        lats = lats[~np.isnan(lats)]

        clat = 0.5 * (np.nanmax(lats) + np.nanmin(lats))
        clon = 0.5 * (np.nanmax(lons) + np.nanmin(lons))
    else:
        # For a point source, it is just the epicenter
        clat = origin.lat
        clon = origin.lon
    return (clon, clat)


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
    Determine if point is located in the US stable tectonic region. This uses
    STREC but only makes use of the stable tectonic region. Any location that
    is not mapped as stable is classified as active.

    Args:
        lon (float): Lognitude.
        lat (float): Latitude.
    Returns:
        bool: Is the point classified as tectonically stable.
    """
    reg = Regionalizer.load()
    region_info = reg.getRegions(lat, lon, 0)
    return region_info['TectonicRegion'] == 'Stable'
