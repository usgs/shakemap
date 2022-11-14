import numpy as np
import logging

from openquake.hazardlib import imt
from openquake.hazardlib import const
from openquake.hazardlib.geo.utils import OrthographicProjection
from openquake.hazardlib.gsim.base import SitesContext
from openquake.hazardlib.gsim.base import DistancesContext
from openquake.hazardlib.gsim.base import RuptureContext
from strec.gmreg import Regionalizer
from impactutils.rupture.edge_rupture import EdgeRupture
from impactutils.rupture.quad_rupture import QuadRupture
from impactutils.rupture.base import Rupture

from shakelib.multigmpe import set_sites_depth_parameters
from shakelib.station import StationList


DEFAULT_ACTIVE_COEFFS = [27.24, 250.4, 579.1]
DEFAULT_STABLE_COEFFS = [63.4, 465.4, 581.3]


def replace_dyfi(stationfile, dyfi_xml, min_nresp=3):
    """Remove any non-instrumented data from station file, add DYFI.

    Args:
        stationfile (str): Existing station data file, presumed to
                           contain old DYFI data.
        dyfi_xml (str): DYFI XML data file, which will be added to
                        instrumented station data.
        min_nresp (int): The minimum number of DYFI responses required for an
                         observation to be added to the stationlist. Default = 3.
    Returns:
        StationList: Object containing merged data.

    """
    stations = StationList.loadFromFiles([stationfile], min_nresp=1)
    # reach into the internal database and find the instrumented stations
    conn = stations.db
    cursor = stations.cursor
    query1 = (
        "SELECT id from station WHERE instrumented = "
        '0 and (network = "DYFI" or network="CIIM")'
    )
    cursor.execute(query1)
    rows = cursor.fetchall()
    for row in rows:
        sid = row[0]
        query2 = f'DELETE FROM amp WHERE station_id="{sid}"'
        cursor.execute(query2)
        conn.commit()
        query3 = f'DELETE FROM station where id="{sid}"'
        cursor.execute(query3)
        conn.commit()

    # now insert the dyfi data
    stations.addData([dyfi_xml], min_nresp)
    return stations


def get_extent(config, ipe, rupture=None):
    """
    Method to compute map extent from rupture. There are numerous methods for
    getting the extent:
        - It can be specified directly in the config file,
        - it can be hard coded for specific magnitude ranges in the config
          file, or
        - it can be based on the MultiGMPE for the event.

    All methods except for the first require a rupture object.

    If no config is provided then a rupture is required and the extent is based
    on a generic set of active/stable.

    Args:
        config (ConfigObj):
            ShakeMap config object.
        ipe (VirtualIPE):
            An VirtualIPE instance.
        rupture (Rupture):
            A ShakeMap Rupture instance.

    Returns:
        tuple: lonmin, lonmax, latmin, latmax rounded outward to the nearest
        30 arc seconds.

    """

    # -------------------------------------------------------------------------
    # Check to see what parameters are specified in the extent config
    # -------------------------------------------------------------------------
    spans = {}
    bounds = []
    offsets = None
    if "extent" in config:
        if "magnitude_spans" in config["extent"]:
            if len(config["extent"]["magnitude_spans"]):
                if isinstance(config["extent"]["magnitude_spans"], dict):
                    spans = config["extent"]["magnitude_spans"]
        if "bounds" in config["extent"]:
            if "extent" in config["extent"]["bounds"]:
                if config["extent"]["bounds"]["extent"][0] != -999.0:
                    bounds = config["extent"]["bounds"]["extent"]
        if "relative_offset" in config["extent"]:
            if isinstance(config["extent"]["relative_offset"], list):
                offsets = config["extent"]["relative_offset"]

    # -------------------------------------------------------------------------
    # Simplest option: extent was specified in the config, use that and exit.
    # -------------------------------------------------------------------------
    if len(bounds):
        xmin, ymin, xmax, ymax = bounds
        return (
            thirty_sec_min(xmin),
            thirty_sec_max(xmax),
            thirty_sec_min(ymin),
            thirty_sec_max(ymax),
        )

    if not rupture or not isinstance(rupture, Rupture):
        raise TypeError(
            "get_extent() requires a rupture object if the extent "
            "is not specified in the config object."
        )

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
        extent = _get_extent_from_multigmpe(rupture, config, ipe)

    if offsets is None:
        xmin, xmax, ymin, ymax = extent
        return (
            thirty_sec_min(xmin),
            thirty_sec_max(xmax),
            thirty_sec_min(ymin),
            thirty_sec_max(ymax),
        )

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

    return (
        thirty_sec_min(xmin),
        thirty_sec_max(xmax),
        thirty_sec_min(ymin),
        thirty_sec_max(ymax),
    )


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
            ymin = clat - span[2] / 2
            ymax = clat + span[2] / 2
            xmin = clon - span[3] / 2
            xmax = clon + span[3] / 2
            break
    if xmin is not None:
        return (xmin, xmax, ymin, ymax)
    return None


def _get_extent_from_multigmpe(rupture, config, ipe):
    """
    Use MultiGMPE to determine extent

    Args:
        rupture (Rupture):
            A ShakeMap Rupture instance.
        config (ConfigObj):
            ShakeMap config object.
        ipe (VirtualIPE):
            An VirtualIPE instance.
    """
    (clon, clat) = _rupture_center(rupture)
    origin = rupture.getOrigin()

    min_mmi = config["extent"]["mmi"]["threshold"]
    sd_types = [const.StdDev.TOTAL]

    # Distance context
    dx = DistancesContext()
    size = 2000
    # This imposes minimum/ maximum distances of:
    #   80 and 800 km; could make this configurable
    d_min = config["extent"]["mmi"]["mindist"]
    d_max = config["extent"]["mmi"]["maxdist"]
    dx.rjb = np.logspace(np.log10(d_min), np.log10(d_max), size)
    dx.rrup = np.sqrt(dx.rjb ** 2 + origin.depth ** 2)
    dx.rhypo = dx.rrup
    dx.repi = dx.rjb
    dx.rx = np.zeros_like(dx.rjb)
    dx.ry0 = np.zeros_like(dx.rjb)
    dx.rvolc = np.zeros_like(dx.rjb)

    # Sites context
    sx = SitesContext()
    # Set to soft soil conditions
    sx.sids = np.array(range(size))
    sx.vs30 = np.full(size, 180.0)
    set_sites_depth_parameters(sx, ipe)
    sx.vs30measured = np.full(size, False, dtype=bool)
    sx.backarc = np.full(size, False, dtype=bool)

    # Rupture context
    rx = RuptureContext()
    rx.mag = origin.mag
    rx.rake = 0.0
    # From WC94...
    # rx.width = 10 ** (-0.76 + 0.27 * rx.mag)
    # rx.dip = 90.0
    # rx.ztor = origin.depth
    # Using parameters from the point rupture...
    rx.width = rupture.getWidth()
    rx.dip = rupture.getDip()
    rx.ztor = rupture.getDepthToTop()
    rx.hypo_depth = origin.depth

    mmi = imt.from_string("MMI")
    imt_mean, _ = ipe.get_mean_and_stddevs(sx, rx, dx, mmi, sd_types)

    # Minimum distance that exceeds threshold MMI?
    dists_exceed_mmi = dx.rjb[imt_mean > min_mmi]
    if len(dists_exceed_mmi):
        mindist_km = np.max(dists_exceed_mmi)
    else:
        mindist_km = d_min

    # Get a projection
    proj = OrthographicProjection(clon - 4, clon + 4, clat + 4, clat - 4)
    if isinstance(rupture, (QuadRupture, EdgeRupture)):
        ruptx, rupty = proj(
            rupture.lons[~np.isnan(rupture.lons)], rupture.lats[~np.isnan(rupture.lats)]
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
    logging.debug(f"Extent: {lonmin[0]:f}, {lonmax[0]:f}, {latmin[0]:f}, {latmax[0]:f}")
    return (
        _round_coord(lonmin[0]),
        _round_coord(lonmax[0]),
        _round_coord(latmin[0]),
        _round_coord(latmax[0]),
    )


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


def thirty_sec_min(coord):
    """
    Round a number to the floor of 30 arc-seconds
    """
    return np.floor(coord * 120.0) / 120.0


def thirty_sec_max(coord):
    """
    Round a number to the ceiling of 30 arc-seconds
    """
    return np.ceil(coord * 120.0) / 120.0


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
    return region_info["TectonicRegion"] == "Stable"
