#!/usr/bin/env python

# stdlib modules
import os
import copy
import warnings
from abc import ABC
from abc import abstractmethod
from abc import abstractproperty
import json
import itertools as it
import re

# third party imports
import numpy as np
import pandas as pd
import scipy.interpolate as spint
from openquake.hazardlib.geo.mesh import Mesh
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.geo.utils import get_orthographic_projection
from openquake.hazardlib.geo import geodetic
from openquake.hazardlib.gsim import base

from shakemap.utils.ecef import latlon2ecef
from shakemap.utils.ecef import ecef2latlon
from shakemap.utils.vector import Vector
from shakemap.utils.timeutils import ShakeDateTime
from shakemap.utils.exception import ShakeMapException
from shakemap.plotting.plotrupture import plot_rupture_wire3d
from shakemap.plotting.plotrupture import map_rupture



#-------------------------------------------------------------------------------
# CONSTANTS and DEFAULTS

# Depth tolerance in km (for determining if top and bottom edges are horizontal)
DEPTH_TOL = 0.05

# Maximum ratio of distance off of the plane (relative to edge length) for the 
# 4th point to be before being considered non-co-planar and adjusted to actually
# be on the plane?
OFFPLANE_TOLERANCE = 0.05

RAKEDICT = {'SS': 0.0, 'NM': -90.0, 'RS': 90.0, 'ALL': None}

DEFAULT_MECH = 'ALL'
DEFAULT_STRIKE = 0.0
DEFAULT_DIP = 90.0
DEFAULT_RAKE = 0.0
DEFAULT_WIDTH = 0.0
DEFAULT_ZTOR = 0.0

#-------------------------------------------------------------------------------




def read_rupture_file(origin, file = None, mesh_dx = 0.5):
    """
    This is a module-level function to read in a rupture file. This allows for
    the ShakeMap 3 text file specification or the ShakeMap 4 JSON rupture format.
    The ShakeMap 3 (".txt" extension) only supports QuadRupture style rupture
    representation and so this method will always return a QuadRupture instance. 
    The ShakeMap 4 JSON format supports QuadRupture and EdgeRupture
    represenations and so this method detects the rupture class and returns the
    appropriate Rupture subclass instance.

    If file is None (default) then it returns a PointRupture.

    Args:
        origin (Origin): A ShakeMap origin instance; required because 
            hypocentral/epicentral distances are computed from the Rupture
            class.
        file (srt): Path to rupture file (optional).
        mesh_dx (float): Target spacing (in km) for rupture discretization;
            default is 0.5 km and it is only used if the rupture file is an 
            EdgeRupture.

    Returns:
        Rupture subclass instance. 

    """
    if file is not None:
        try:
            #-----------------------------------------------------------------------
            # First, try to read as a json file
            #-----------------------------------------------------------------------
            if isinstance(file, str):
                with open(file) as f:
                    d = json.load(f)
            else:
                d = json.loads(str(file))

            rupt = json_to_rupture(d, origin, mesh_dx = mesh_dx)

        except json.JSONDecodeError:
            #-----------------------------------------------------------------------
            # Reading as json failed, so hopefully it is a ShakeMap 3 text file
            #-----------------------------------------------------------------------
            try:
                d = text_to_json(file)
                rupt = json_to_rupture(d, origin, mesh_dx = mesh_dx)
            except:
                raise Exception("Unknown rupture file format.")
    else:
        if origin is None:
            raise Exception("Origin requred if no rupture file is provided.")
        rupt = PointRupture(origin)
    return rupt




def json_to_rupture(d, origin, mesh_dx = 0.5):
    """
    Method returns either a QuadRupture or EdgeRupture object based on a 
    GeoJSON dictionary. 

    Args: 
        d (dict): Rupture GeoJSON dictionary.
        origin (Origin): A ShakeMap origin object.
        mesh_dx (float): Target spacing (in km) for rupture discretization;
            default is 0.5 km and it is only used if the rupture file is an
            EdgeRupture.

    Returns:
        a Rupture subclass.

    """
    validate_json(d)

    # Is this a QuadRupture or an EdgeRupture?
    valid_quads = is_quadrupture_class(d)

    if valid_quads is True:
        rupt = QuadRupture(d, origin)
    else:
        rupt = EdgeRupture(d, origin, mesh_dx = mesh_dx)

    return rupt


def text_to_json(file):
    """
    Read in old ShakeMap 3 textfile rupture format and convert to GeoJSON. 

    Args:
        rupturefile (srt): Path to rupture file OR file-like object in GMT
            psxy format, where:

                * Rupture vertices are space separated lat, lon, depth triplets 
                  on a single line.
                * Rupture groups are separated by lines containing ">"
                * Rupture groups must be closed.
                * Verticies within a rupture group must start along the top edge
                  and move in the strike direction then move to the bottom edge
                  and move back in the opposite direction.

    Returns:
        dict: GeoJSON rupture dictionary.

    """

    #---------------------------------------------------------------------------
    # First read in the data
    #---------------------------------------------------------------------------
    x = []
    y = []
    z = []
    isFile = False
    if isinstance(file, str):
        isFile = True
        file = open(file, 'rt')
        lines = file.readlines()
    else:
        lines = file.readlines()
    reference = ''
    for line in lines:
        sline = line.strip()
        if sline.startswith('#'):
            reference += sline
            continue
        if sline.startswith('>'):
            if len(x):  # start of new line segment
                x.append(np.nan)
                y.append(np.nan)
                z.append(np.nan)
                continue
            else:  # start of file
                continue
        if not len(sline.strip()):
            continue
        parts = sline.split()
        if len(parts) < 3:
            raise ShakeMapException(
                'Rupture file %s has no depth values.' % file)
        y.append(float(parts[0]))
        x.append(float(parts[1]))
        z.append(float(parts[2]))
    if isFile:
        file.close()

    # Construct GeoJSON dictionary

    coords = []
    poly = []
    for lon, lat, dep in zip(x, y, z):
        if np.isnan(lon):
            coords.append(poly)
            poly = []
        else:
            poly.append([lon, lat, dep])
    if poly != []:
        coords.append(poly)

    d = {
        "type": "FeatureCollection",
        "metadata": {},
        "features":[
            {
                "type": "Feature",
                "properties": {
                    "rupture type": "rupture extent",
		    "reference": reference
                },
                "geometry": {
	            "type": "MultiPolygon",
	            "coordinates":[coords]
                }
            }
        ]
    }
    return d



def validate_json(d):
    """
    Check that the JSON format is acceptable. This is only for requirements that
    are common to both QuadRupture and EdgeRupture.

    Args:
        d (dict): Rupture JSON dictionary.
    """
    if d['type'] != 'FeatureCollection':
        raise Exception('JSON file is not a \"FeatureColleciton\".')

    if len(d['features']) != 1:
        raise Exception('JSON file should contain excactly one feature.')

    f = d['features'][0]

    if 'reference' not in f['properties'].keys():
        raise Exception('Feature property dictionary should contain '\
                        '\"referencey\" key.')

    if f['type'] != 'Feature':
        raise Exception('Feature type should be \"Feature\".')

    geom = f['geometry']

    if geom['type'] != 'MultiPolygon':
        raise Exception('Geometry type should be \"MultiPolygon\".')

    if 'coordinates' not in geom.keys():
        raise Exception('Geometry dictionary should contain \"coordinates\" '\
                        'key.')

    polygons = geom['coordinates'][0]

    n_polygons = len(polygons)
    for i in range(n_polygons):
        p = polygons[i]
        n_points = len(p)
        if n_points % 2 == 0:
            raise Exception('Number of points in polyon must be odd.')

        if p[0] != p[-1]:
            raise Exception('First and last points in polygon must be '\
                            'identical.')

        n_pairs = int((n_points - 1)/2)
        for j in range(n_pairs):
            #-------------------------------------------------------------------
            # Points are paired and in each pair the top is first, as in:
            #
            #      _.-P1-._
            #   P0'        'P2---P3
            #   |                  \
            #   P7---P6----P5-------P4
            #
            # Pairs: P0-P7, P1-P6, P2-P5, P3-P4
            #-------------------------------------------------------------------
            top_depth = p[j][2]
            bot_depth = p[-(j+2)][2]
            if top_depth > bot_depth:
                raise Exception('Top points must be ordered before bottom points.')


def is_quadrupture_class(d):
    """
    Check if JSON file fulfills QuadRupture class criteria:
    
        - Are top and bottom edges horizontal?
        - Are the four points in each quad coplanar?

    Args:
        d (dict): Rupture JSON dictionary.

    Returns:
        bool: Can the rupture be represented in the QuadRupture class?
    """
    isQuad = True

    f = d['features'][0]
    geom = f['geometry']
    polygons = geom['coordinates'][0]
    n_polygons = len(polygons)
    for i in range(n_polygons):
        p = polygons[i]
        n_points = len(p)
        n_pairs = int((n_points - 1)/2)

        # Within each polygon, top and bottom edges must be horizontal
        depths = [pt[2] for pt in p]
        tops = np.array(depths[0:n_pairs])
        if not np.isclose(tops[0], tops, rtol = 0, atol = DEPTH_TOL).all():
            isQuad = False
        bots = np.array(depths[(n_pairs):-1])
        if not np.isclose(bots[0], bots, rtol = 0, atol = DEPTH_TOL).all():
            isQuad = False

        n_quads = n_pairs - 1
        for j in range(n_quads):
            # Four points of each quad should be co-planar within a tolerance
            quad = [Point(p[j][0], p[j][1], p[j][2]),
                    Point(p[j+1][0], p[j+1][1], p[j+1][2]),
                    Point(p[-(j+3)][0], p[-(j+3)][1], p[-(j+3)][2]),
                    Point(p[-(j+2)][0], p[-(j+2)][1], p[-(j+2)][2])]

            test = is_quad(quad)
            if test[0] is False:
                isQuad = False

    return isQuad


def is_quad(q):
    """
    Checks that an individual quad is coplanar. 

    Args: 
        q (list): A quadrilateral; list of four OQ Points.

    Returns:
        tuple: Bool for whether or not the points are planar within tolerance;
            and also the corrected quad where p2 is adjusted to be on the same
            plane as the other points.
    """
    P0, P1, P2, P3 = q

    # Convert points to ECEF
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)

    # Unit vector along top edge
    v0 = (p1 - p0).norm()

    # Distance along bottom edge
    d = (p3 - p2).mag()

    # New location for p2 by extending from p3 the same distance and
    # direction that p1 is from p0:
    new_p2 = p3 + v0*d

    # How far off of the plane is the origin p2?
    planepoints = [p0, p1, p2]
    dist = get_distance_to_plane(planepoints, p2)

    # Is it close enough?
    if dist / d > OFFPLANE_TOLERANCE:
        on_plane = False
    else:
        on_plane = True

    # Fixed quad
    fquad = [p0.toPoint(),
             p1.toPoint(),
             new_p2.toPoint(),
             p3.toPoint()]

    return (on_plane, fquad)



class Rupture(ABC):
    """
    Abstract base class for ruptures.

    Note on terminology:

        - There are three Ruptuer subclasses: PointRupture, QuadRupture, and
          EdgeRupture.
        - PointRupture represents the rupture as a point source.
        - QuadRupture and EdgeRupture are two different finite source 
          representations.
        - A finite rupture is composed of segments. For QuadRupture, a segment
          is a quadrilaterial; for an EdgeRupture, a segment is a line
          connecting two points.
        - Segments are grouped with a common "group index".
        - Segments within a group must be continuous.
        - The QuadRupture class requires that each segment is a quadrilateral
          with horizontal top and obttom edges.
        - The EdgeRupture class allows for arbitrarily complex top and bottom
          edge specification. 

    """


    def writeGeoJson(self, file):
        """
        Write the rupture to a GeoJson file.

        Args:
            file (str): Name of file. 
        """
        with open(file, 'w') as f:
            json.dump(self._geojson, f)


    @abstractmethod
    def getLength(self):
        """
        Returns:
            float: Rupture length in km.
        """
        pass

    @abstractmethod
    def getWidth(self):
        """
        Returns:
            float: Rupture width in km.
        """
        pass

    @abstractmethod
    def getArea(self):
        """
        Returns:
            float: Rupture area in square km.
        """
        pass

    @abstractmethod
    def getStrike(self):
        """
        Return strike angle. If rupture consists of multiple quadrilaterals, the
        average strike angle, weighted by quad length, is returned.
        Note: for ruptures with quads where the strike angle changes by 180 deg
        due to reverses in dip direction are problematic and not handeled well
        by this algorithm.

        Returns:
            float: Strike angle in degrees.
        
        """
        pass

    @abstractmethod
    def getDip(self):
        pass

    @abstractmethod
    def getDepthToTop(self):
        """
        Returns:
           float: Average dip in degrees.

        """
        pass

    @abstractmethod
    def getQuadrilaterals(self):
        """
        Method to return rupture quadrilaterals. Returns None for
        PointRupture.
        """
        pass

    def getReference(self):
        """
        Returns:
           float: Average dip in degrees.

        """
        return self._reference

    @abstractproperty
    def lats(self):
        pass

    @abstractproperty
    def lons(self):
        pass

    @abstractproperty
    def depths(self):
        pass

    def getRuptureContext(self, gmpelist):
        """
        Returns:
            An Openquake 
        `RuptureContext <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#openquake.hazardlib.gsim.base.RuptureContext>`__.

        Args:
            gmpelist (list): List of hazardlib GMPE objects.

        Returns:
            RuptureContext object with all known parameters filled in.

        """

        origin = self._origin

        # rupturecontext constructor inputs:
        # 'mag', 'strike', 'dip', 'rake', 'ztor', 'hypo_lon', 'hypo_lat',
        # 'hypo_depth', 'width', 'hypo_loc'

        rx = base.RuptureContext()
        rx.mag = origin.mag
        rx.strike = self.getStrike()
        rx.dip = self.getDip()
        rx.ztor = self.getDepthToTop()
        rx.width = self.getWidth()

        if hasattr(origin, 'rake'):
            rx.rake = origin.rake
        elif hasattr(origin, 'mech'):
            mech = origin.mech
            rx.rake = RAKEDICT[mech]
        else:
            rx.rake = RAKEDICT['ALL']

        rx.hypo_lat = origin.lat
        rx.hypo_lon = origin.lon
        rx.hypo_depth = origin.depth

        return rx

    def computeRhyp(self, lon, lat, depth):
        """
        Method for computing hypocentral distance. 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Hypocentral distance (km).
        """
        origin = self._origin
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        rhyp = geodetic.distance(origin.lon, origin.lat, origin.depth,
                                 lon, lat, depth)
        rhyp = rhyp.reshape(oldshape)
        return rhyp

    def computeRepi(self, lon, lat, depth):
        """
        Method for computing epicentral distance. 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Epicentral distance (km).
        """
        origin = self._origin
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        repi = geodetic.distance(origin.lon, origin.lat, 0.0,
                                 lon, lat, depth)
        repi = repi.reshape(oldshape)
        return repi        


    @abstractmethod
    def computeRjb(self, lon, lat, depth):
        """
        Method for computing Joyner-Boore distance. 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Joyner-Boore distance (km).

        """
        pass

    @abstractmethod
    def computeRrup(self, lon, lat, depth):
        """
        Method for computing rupture distance. 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Rupture distance (km).

        """
        pass

    @abstractmethod
    def computeGC2(self, lon, lat, depth):
        """
        Method for computing version 2 of the Generalized Coordinate system 
        (GC2) by Spudich and Chiou OFR 2015-1028. 
        
        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            dict: Dictionary with keys for each of the GC2-related distances,
                which include 'rx', 'ry', 'ry0', 'U', 'T'. 
        """
        pass

    def plot3d(self):
        """
        Method for making a quick 3D wireframe plot of rupture.
        """
        plot_rupture_wire3d(self)

    def map(self):
        """
        Method for making a quick map of the fault.
        """
        map_rupture(self)


class PointRupture(Rupture):
    """
    Rupture class for point sources. The purpose is to gracefully handle:

        - Requests for rupture distances when no rupture is available. 
        - Provide reasonable default values for rupture parameters. 
    """

    def __init__(self, origin, reference = ""):
        """
        Constructs a PointRupture instance. 

        Args:
            origin (Origin): Reference to a ShakeMap Origin instance.
            reference (str): Citable reference for rupture. 

        Returns:
            PointRupture instance.
        """
        self._origin = origin
        self._reference = reference

        coords = [origin.lon, origin.lat, origin.depth]

        d = {"type":"FeatureCollection", 
             "metadata":{},
             "features":[{
                    "type":"Feature",
                    "properties":{
                        "rupture type":"rupture extent",
                        "reference":reference,
                        },
                    "geometry":{
                        "type":"Point",
                        "coordinates":coords
                        }
            }]}


        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, ShakeDateTime):
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                d['metadata'][k] = v

        self._geojson = d


    def getLength(self):
        """
        Rupture length, which is None for a PointRupture.
        Could potentially put in a default value based on magnitude, and add
        an Origin argument.
        """
        return None

    def getWidth(self):
        """
        Rupture width.
        Could potentially put in a default value based on magnitude, and add
        an Origin argument.
        """
        return DEFAULT_WIDTH

    def getArea(self):
        """
        Rupture area, which is None for a PointRupture.
        Could potentially put in a default value based on magnitude, and add
        an Origin argument.
        """
        return None

    def getStrike(self):
        """
        Strike, which is None. 
        Could potentially get from strec or something?
        """
        return DEFAULT_STRIKE

    def getDip(self):
        """
        Dip, which is None.
        Could potentially get from strec or something?
        """
        return DEFAULT_DIP

    def getDepthToTop(self):
        """
        Depth to top of rupture.
        Could get from hypo/magnitude?
        """
        return DEFAULT_ZTOR

    def getQuadrilaterals(self):
        return None

    @property
    def lats(self):
        """
        Returns rupture latitudes, which is just the hypocenter for a PointRupture.
        """
        return self._origin.lat

    @property
    def lons(self):
        """
        Returns rupture longitudes, which is just the hypocenter for a PointRupture.
        """
        return self._origin.lon

    @property
    def depths(self):
        """
        Returns rupture depths, which is just the hypocenter for a PointRupture.
        """
        return self._origin.depth



    def computeRjb(self, lon, lat, depth, var = False):
        """
        Convert point-distances to Joyner-Boore distances based on magnitude.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).
            var (bool): Also return variance of prediction. 

        Returns:
            If var is True then this returns a tuple of two arrays: first, the 
                predicted Rjb values, and second an array of the variance of
                those predictions. If var is False then just the first element
                of the tuple is returned.
        """
        cdir, tmp = os.path.split(__file__)
        origin = self._origin

        # -------------------
        # Sort out file names
        # -------------------
        mech = origin.mech
        if not hasattr(origin, '_tectonic_region'):
            rf = os.path.join(
                cdir, "data", "ps2ff",
                "Rjb_WC94_mechA_ar1p0_seis0_20_Ratios.csv")
            vf = os.path.join(
                cdir, "data", "ps2ff",
                "Rjb_WC94_mechA_ar1p0_seis0_20_Var.csv")
        elif origin._tectonic_region == 'Active Shallow Crust':
            if mech == 'ALL':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_WC94_mechA_ar1p7_seis0_20_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_WC94_mechA_ar1p7_seis0_20_Var.csv")
            elif mech == 'RS':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_WC94_mechR_ar1p7_seis0_20_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_WC94_mechR_ar1p7_seis0_20_Var.csv")
            elif mech == 'NM':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_WC94_mechN_ar1p7_seis0_20_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_WC94_mechN_ar1p7_seis0_20_Var.csv")
            elif mech == 'SS':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_WC94_mechSS_ar1p7_seis0_20_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_WC94_mechSS_ar1p7_seis0_20_Var.csv")
        elif origin._tectonic_region == 'Stable Shallow Crust':
            if mech == 'ALL':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_S14_mechA_ar1p0_seis0_15_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_S14_mechA_ar1p0_seis0_15_Var.csv")
            elif mech == 'RS':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_S14_mechR_ar1p0_seis0_15_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_S14_mechR_ar1p0_seis0_15_Var.csv")
            elif mech == 'NM':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_S14_mechN_ar1p0_seis0_15_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_S14_mechN_ar1p0_seis0_15_Var.csv")
            elif mech == 'SS':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_S14_mechSS_ar1p0_seis0_15_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rjb_S14_mechSS_ar1p0_seis0_15_Var.csv")
        else:
            warnings.warn(
                'Unsupported tectonic region; using coefficients for unknown'
                'tectonic region.')
            rf = os.path.join(
                cdir, "data", "ps2ff",
                "Rjb_WC94_mechA_ar1p0_seis0_20_Ratios.csv")
            vf = os.path.join(
                cdir, "data", "ps2ff",
                "Rjb_WC94_mechA_ar1p0_seis0_20_Var.csv")

        # -----------------
        # Start with ratios
        # -----------------
        repi2rjb_ratios_tbl = pd.read_csv(rf, comment='#')
        r2rrt_cols = repi2rjb_ratios_tbl.columns[1:]
        mag_list = []
        for column in (r2rrt_cols):
            if re.search('R\d+\.*\d*', column):
                magnitude = float(re.findall(
                    'R(\d+\.*\d*)', column)[0])
                mag_list.append(magnitude)
        mag_list = np.array(mag_list)
        dist_list = np.log(np.array(repi2rjb_ratios_tbl['Repi_km']))
        repi2rjb_grid = repi2rjb_ratios_tbl.values[:, 1:]
        repi2rjb_obj = spint.RectBivariateSpline(
            dist_list, mag_list, repi2rjb_grid, kx=1, ky=1)

        def repi2rjb_tbl(repi, M):
            ratio = repi2rjb_obj.ev(np.log(repi), M)
            rjb = repi * ratio
            return rjb

        repis = self.computeRepi(lon, lat, depth)
        mags = np.ones_like(repis) * origin.mag
        rjb_hat = repi2rjb_tbl(repis, mags)

        # -------------------
        # Additional Variance
        # -------------------
        repi2rjbvar_ratios_tbl = pd.read_csv(vf, comment='#')
        repi2rjbvar_grid = repi2rjbvar_ratios_tbl.values[:, 1:]
        repi2rjbvar_obj = spint.RectBivariateSpline(
            dist_list, mag_list, repi2rjbvar_grid, kx=1, ky=1)
        rjbvar = repi2rjbvar_obj.ev(np.log(repis), mags)

        if var is True:
            return (rjb_hat, rjbvar)
        else:
            return rjb_hat




    def computeRrup(self, lon, lat, depth, var = False):
        """
        Convert point-distances to rupture distances based on magnitude.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).
            var (bool): Also return variance of prediction.

        Returns:
            If var is True then this returns a tuple of two arrays: first, the
                predicted Rjb values, and second an array of the variance of
                those predictions. If var is False then just the first element
                of the tuple is returned.
        """
        cdir, tmp = os.path.split(__file__)
        origin = self._origin

        # -------------------
        # Sort out file names
        # -------------------
        rake = float(origin.rake)
        mech = rake_to_mech(rake)
        if not hasattr(origin, '_tectonic_region'):
            rf = os.path.join(
                cdir, "data", "ps2ff",
                "Rrup_WC94_mechA_ar1p0_seis0-20_Ratios.csv")
            vf = os.path.join(
                cdir, "data", "ps2ff",
                "Rrup_WC94_mechA_ar1p0_seis0-20_Var.csv")
        elif origin._tectonic_region == 'Active Shallow Crust':
            if mech == 'ALL':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_WC94_mechA_ar1p7_seis0-20_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_WC94_mechA_ar1p7_seis0-20_Var.csv")
            elif mech == 'RS':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_WC94_mechR_ar1p7_seis0-20_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_WC94_mechR_ar1p7_seis0-20_Var.csv")
            elif mech == 'NM':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_WC94_mechN_ar1p7_seis0-20_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_WC94_mechN_ar1p7_seis0-20_Var.csv")
            elif mech == 'SS':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_WC94_mechSS_ar1p7_seis0-20_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_WC94_mechSS_ar1p7_seis0-20_Var.csv")
        elif origin._tectonic_region == 'Stable Shallow Crust':
            if mech == 'ALL':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_S14_mechA_ar1p0_seis0-15_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_S14_mechA_ar1p0_seis0-15_Var.csv")
            elif mech == 'RS':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_S14_mechR_ar1p0_seis0-15_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_S14_mechR_ar1p0_seis0-15_Var.csv")
            elif mech == 'NM':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_S14_mechN_ar1p0_seis0-15_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_S14_mechN_ar1p0_seis0-15_Var.csv")
            elif mech == 'SS':
                rf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_S14_mechSS_ar1p0_seis0-15_Ratios.csv")
                vf = os.path.join(
                    cdir, "data", "ps2ff",
                    "Rrup_S14_mechSS_ar1p0_seis0-15_Var.csv")
        else:
            warnings.warn(
                'Unsupported tectonic region; using coefficients for unknown'
                'tectonic region.')
            rf = os.path.join(
                cdir, "data", "ps2ff",
                "Rrup_WC94_mechA_ar1p0_seis0-20_Ratios.csv")
            vf = os.path.join(
                cdir, "data", "ps2ff",
                "Rrup_WC94_mechA_ar1p0_seis0-20_Var.csv")

        # -----------------
        # Start with ratios
        # -----------------
        repi2rrup_ratios_tbl = pd.read_csv(rf, comment='#')
        r2rrt_cols = repi2rrup_ratios_tbl.columns[1:]
        mag_list = []
        for column in (r2rrt_cols):
            if re.search('R\d+\.*\d*', column):
                magnitude = float(re.findall(
                    'R(\d+\.*\d*)', column)[0])
                mag_list.append(magnitude)
        mag_list = np.array(mag_list)
        dist_list = np.log(np.array(repi2rrup_ratios_tbl['Repi_km']))
        repi2rrup_grid = repi2rrup_ratios_tbl.values[:, 1:]
        repi2rrup_obj = spint.RectBivariateSpline(
            dist_list, mag_list, repi2rrup_grid, kx=1, ky=1)

        def repi2rrup_tbl(repi, M):
            ratio = repi2rrup_obj.ev(np.log(repi), M)
            rrup = repi * ratio
            return rrup

        repis = self.computeRepi(lon, lat, depth)
        mags = np.ones_like(repis) * origin.mag
        rrup_hat = repi2rrup_tbl(repis, mags)

        # -------------------
        # Additional Variance
        # -------------------
        repi2rrupvar_ratios_tbl = pd.read_csv(vf, comment='#')
        repi2rrupvar_grid = repi2rrupvar_ratios_tbl.values[:, 1:]
        repi2rrupvar_obj = spint.RectBivariateSpline(
            dist_list, mag_list, repi2rrupvar_grid, kx=1, ky=1)
        rrupvar = repi2rrupvar_obj.ev(np.log(repis), mags)

        if var is True:
            return (rrup_hat, rrupvar)
        else:
            return rrup_hat

    def computeGC2(self, lon, lat, depth):
        """
        Method for computing version 2 of the Generalized Coordinate system 
        (GC2) by Spudich and Chiou OFR 2015-1028. 
        
        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            dict: Dictionary with keys for each of the GC2-related distances,
                which include 'rx', 'ry', 'ry0', 'U', 'T'. 
        """
        # This just returns defaults of zero, which will hopefully behave gracefully
        # as used in GMPEs.
        dict = {"rx":np.zeros_like(lon),
                "ry":np.zeros_like(lon),
                "ry0":np.zeros_like(lon),
                "U":np.zeros_like(lon),
                "T":np.zeros_like(lon) 
                }
        return dict
        

class EdgeRupture(Rupture):
    """
    Rupture class that representst the rupture surface by specifying the top
    edge and the bottom edges. These edges do not need to be horizontal. The
    freedom to allow for non-horizontal edges (as compared to QuadRupture)
    comes at the cost of slower distance calculations. This is because the 
    rupture must be discretized and then the distances are compued in a brute
    force fashion based on this mesh, which can be quite large. 
    """

    def __init__(self, d, origin, mesh_dx = 0.5):
        """
        Initialization of an EdgeRupture from a GeoJSON dictionary and an Origin.

        Args: 
            d (dict): Rupture GeoJSON dictionary.
            origin (Origin): Reference to a ShakeMap Origin object.
            mesh_dx (float): Target spacing (in km) for rupture discretization;
                default is 0.5 km and it is only used if the rupture file is an
                EdgeRupture.

        Returns: 
            EdgeRupture instance.

        """
        polys = d['features'][0]['geometry']['coordinates'][0]
        n_polygons = len(polys)
        toplons = np.empty(shape=(0, 0))
        toplats = np.empty(shape=(0, 0))
        topdeps = np.empty(shape=(0, 0))
        botlons = np.empty(shape=(0, 0))
        botlats = np.empty(shape=(0, 0))
        botdeps = np.empty(shape=(0, 0))
        g_ind = 0
        group_index = []
        for i in range(n_polygons):
            p = polys[i]
            n_points = len(p)
            n_pairs = int((n_points - 1)/2)
            
            p_lons = [pt[0] for pt in p][0:-1]
            p_lats = [pt[1] for pt in p][0:-1]
            p_depths = [pt[2] for pt in p][0:-1]

            tlon = np.array(p_lons[0:n_pairs])
            blon = np.array(p_lons[(n_pairs):])[::-1]
            tlat = np.array(p_lats[0:n_pairs])
            blat = np.array(p_lats[(n_pairs):])[::-1]
            tdep = np.array(p_depths[0:n_pairs])
            bdep = np.array(p_depths[(n_pairs):])[::-1]
            
            toplons = np.append(toplons, tlon)
            toplats = np.append(toplats, tlat)
            topdeps = np.append(topdeps, tdep)
            botlons = np.append(botlons, blon)
            botlats = np.append(botlats, blat)
            botdeps = np.append(botdeps, bdep)

            group_index.extend([g_ind]*n_pairs)
            g_ind = g_ind + 1

        reference = d['features'][0]['properties']['reference']

        self._geojson = d

        self._toplons = np.array(toplons)
        self._toplats = np.array(toplats)
        self._topdeps = np.array(topdeps)
        self._botlons = np.array(botlons)
        self._botlats = np.array(botlats)
        self._botdeps = np.array(botdeps)
        self._origin = origin
        self._group_index = np.array(group_index)
        self._mesh_dx = mesh_dx
        self._reference = reference
        self._computeStikeDip()


    @classmethod
    def fromArrays(cls, toplons, toplats, topdeps, botlons, botlats, botdeps,
                   origin, group_index = None, mesh_dx = 0.5, reference = ''):
        """
        Class method to initialize an EdgeRupture class from arrays of 
        longitude, latitude, and depth along the top and bottom edges. 

        Args:
            toplons (ndarray): Array of top edge longitudes.
            toplats (ndarray): Array of top edge latitudes. 
            topdeps (ndarray): Array of top edge depths (km).
            botlons (ndarray): Array of bot edge longitudes.
            botlats (ndarray): Array of bot edge latitudes. 
            botdeps (ndarray): Array of bot edge depths (km).
            origin (Origin): Reference to a ShakeMap Origin object. 
            group_index (ndarray): Optional array of group index. 
                If None, then assume only single group. 
            mesh_dx (float): Target spacing (in km) for rupture discretization.
            reference (str): Citable reference for rupture.
 
        Returns: 
            EdgeRupture instance.

        """
        toplons = np.array(toplons)
        toplats = np.array(toplats)
        topdeps = np.array(topdeps)
        botlons = np.array(botlons)
        botlats = np.array(botlats)
        botdeps = np.array(botdeps)
        if group_index is not None:
            group_index = np.array(group_index)
        else:
            group_index = np.zeros_like(toplons)

        coords = []
        u_groups = np.unique(group_index)
        n_groups = len(u_groups)
        for i in range(n_groups):
            ind = np.where(u_groups[i] == group_index)[0]
            lons = np.concatenate([toplons[ind], botlons[ind][::-1],
                                   toplons[ind][0].reshape((1,))])
            lats = np.concatenate([toplats[ind], botlats[ind][::-1],
                                   toplats[ind][0].reshape((1,))])
            deps = np.concatenate([topdeps[ind], botdeps[ind][::-1],
                                   topdeps[ind][0].reshape((1,))])
            poly = []
            for lon,lat,dep in zip(lons, lats, deps):
                poly.append([lon, lat, dep])
            coords.append(poly)

        d = {"type":"FeatureCollection", 
             "metadata":{},
             "features":[{
                    "type":"Feature",
                    "properties":{
                        "rupture type":"rupture extent",
                        "reference":reference,
                        },
                    "geometry":{
                        "type":"MultiPolygon",
                        "coordinates":[coords]
                        }
            }]}

        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, ShakeDateTime):
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                d['metadata'][k] = v

        return cls(d, origin, mesh_dx)


    def getLength(self):
        """
        Compute length of rupture (km). For EdgeRupture, we compute the length
        as the length of the top edge projected to the surface.

        Returns:
            float: Rupture length in km.
        """
        lons = self._toplons
        lats = self._toplats
        seg = self._group_index
        groups = np.unique(seg)
        ng = len(groups)
        rlength = 0
        for i in range(ng):
            group_segments = np.where(groups[i] == seg)[0]
            nseg = len(group_segments) - 1
            for j in range(nseg):
                ind = group_segments[j]
                P0 = Point(lons[ind], lats[ind])
                P1 = Point(lons[ind + 1], lats[ind + 1])
                dist = P0.distance(P1)
                rlength = rlength + dist
        return rlength

    def getWidth(self):
        """
        Compute average rupture width in km. For EdgeRupture, we compute this as
        the rupture area divided by teh rupture length. 

        Returns:
            float: Rupture width in km.
        """
        area = self.getArea()
        length = self.getLength()
        return area/length

    def getArea(self):
        """
        Compute the rupture area. For EdgeRupture, we compute this by grouping
        the traces into "quadrilaterals" for which the verticies may not be
        co-planar. We then break up the quadrilaterals into triangles for which
        we can compute area. 

        Returns:
            float: Rupture area in square km.
        """
        seg = self._group_index
        groups = np.unique(seg)
        ng = len(groups)
        area = 0
        for i in range(ng):
            group_segments = np.where(groups[i] == seg)[0]
            nseg = len(group_segments) - 1
            for j in range(nseg):
                ind = group_segments[j]
                p0 = latlon2ecef(self._toplats[ind],
                                 self._toplons[ind],
                                 self._topdeps[ind])
                p1 = latlon2ecef(self._toplats[ind+1],
                                 self._toplons[ind+1],
                                 self._topdeps[ind+1])
                p2 = latlon2ecef(self._botlats[ind+1],
                                 self._botlons[ind+1],
                                 self._botdeps[ind+1])
                p3 = latlon2ecef(self._botlats[ind],
                                 self._botlons[ind],
                                 self._botdeps[ind])
                a = np.sqrt((p1[0] - p0[0])**2 +
                            (p1[1] - p0[1])**2 +
                            (p1[2] - p0[2])**2)
                b = np.sqrt((p2[0] - p0[0])**2 +
                            (p2[1] - p0[1])**2 +
                            (p2[2] - p0[2])**2)
                c = np.sqrt((p2[0] - p1[0])**2 +
                            (p2[1] - p1[1])**2 +
                            (p2[2] - p1[2])**2)
                s = (a + b + c)/2
                A1 = np.sqrt(s*(s - a)*(s - b)*(s - c))
                a = np.sqrt((p0[0] - p3[0])**2 +
                            (p0[1] - p3[1])**2 +
                            (p0[2] - p3[2])**2)
                b = np.sqrt((p2[0] - p3[0])**2 +
                            (p2[1] - p3[1])**2 +
                            (p2[2] - p3[2])**2)
                c = np.sqrt((p0[0] - p2[0])**2 +
                            (p0[1] - p2[1])**2 +
                            (p0[2] - p2[2])**2)
                s = (a + b + c)/2
                A2 = np.sqrt(s*(s - a)*(s - b)*(s - c))
                area = area + (A1 + A2)/1000/1000
        return area


    def getStrike(self):
        """
        Return representative strike for this rupture. Note that strike
        can vary along the rupture.

        Returns:
            float: Strike angle in degrees.
        """
        return self._strike

    def getDip(self):
        """
        Representative dip for this rupture. Note that dip
        can vary along the rupture.

        Returns:
            float: dip angle in degrees.
        """
        return self._dip

    @property
    def lats(self):
        """
        Return an array of latitudes for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected
        segments.

        Returns:
            array: Numpy array of closed-loop latitude values; disconnected
                segments are separated by nans.
        """
        lats = np.empty(shape=(0,))
        groups = self._group_index
        u_groups = np.unique(groups)
        ng = len(u_groups)
        nan = np.array(np.nan).reshape(1,)
        for i in range(ng):
            top_lats = self._toplats[groups == u_groups[i]]
            top0 = top_lats[0].reshape((1,))
            bot_lats = self._botlats[groups == u_groups[i]]
            lats = np.concatenate((lats, top_lats, bot_lats[::-1], top0, nan))
        return np.array(lats)

    @property
    def lons(self):
        """
        Return an array of longitudes for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected
        segments.

        Returns:
            array: Numpy array of closed-loop longitude values; disconnected
                segments are separated by nans.
        """
        lons = np.empty(shape=(0,))
        groups = self._group_index
        u_groups = np.unique(groups)
        ng = len(u_groups)
        nan = np.array(np.nan).reshape(1,)
        for i in range(ng):
            top_lons = self._toplons[groups == u_groups[i]]
            top0 = top_lons[0].reshape((1,))
            bot_lons = self._botlons[groups == u_groups[i]]
            lons = np.concatenate((lons, top_lons, bot_lons[::-1], top0, nan))
        return np.array(lons)

    @property
    def depths(self):
        """
        Return an array of depths for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected
        segments.

        Returns:
            array: Numpy array of closed-loop latitude values; disconnected
                segments are separated by nans.
        """
        deps = np.empty(shape=(0,))
        groups = self._group_index
        u_groups = np.unique(groups)
        ng = len(u_groups)
        nan = np.array(np.nan).reshape(1,)
        for i in range(ng):
            top_deps = self._topdeps[groups == u_groups[i]]
            top0 = top_deps[0].reshape((1,))
            bot_deps = self._botdeps[groups == u_groups[i]]
            deps = np.concatenate((deps, top_deps, bot_deps[::-1], top0, nan))
        return np.array(deps)

    def _getGroupIndex(self):
        """
        Return a list of segment group indexes.

        Returns:
            list: Segment group indexes; length equals the number of quadrilaterals.
        """
        return copy.deepcopy(self._group_index)


    def _computeStikeDip(self):
        """
        Loop over all triangles and get the average normal, north, and up vectors
        in ECEF. Use these to compute a representative strike and dip. 
        """
        seg = self._group_index
        groups = np.unique(seg)
        ng = len(groups)
        norm_vec = Vector(0, 0, 0)
        north_vec = Vector(0, 0, 0)
        up_vec = Vector(0, 0, 0)
        for i in range(ng):
            group_segments = np.where(groups[i] == seg)[0]
            nseg = len(group_segments) - 1
            for j in range(nseg):
                ind = group_segments[j]
                P0 = Point(self._toplons[ind],
                           self._toplats[ind],
                           self._topdeps[ind])
                P1 = Point(self._toplons[ind+1],
                           self._toplats[ind+1],
                           self._topdeps[ind+1])
                P2 = Point(self._botlons[ind+1],
                           self._botlats[ind+1],
                           self._botdeps[ind+1])
                P3 = Point(self._botlons[ind],
                           self._botlats[ind],
                           self._botdeps[ind])
                P1up = Point(self._toplons[ind+1],
                             self._toplats[ind+1],
                             self._topdeps[ind+1]-1.0)
                P1N = Point(self._toplons[ind+1],
                            self._toplats[ind+1]+0.001,
                            self._topdeps[ind+1])
                P3up = Point(self._botlons[ind],
                             self._botlats[ind],
                             self._botdeps[ind]-1.0)
                P3N = Point(self._botlons[ind],
                            self._botlats[ind]+0.001,
                            self._botdeps[ind])
                p0 = Vector.fromPoint(P0)
                p1 = Vector.fromPoint(P1)
                p2 = Vector.fromPoint(P2)
                p3 = Vector.fromPoint(P3)
                p1up = Vector.fromPoint(P1up)
                p1N = Vector.fromPoint(P1N)
                p3up = Vector.fromPoint(P3up)
                p3N = Vector.fromPoint(P3N)

                # Sides
                s01 = p1 - p0
                s02 = p2 - p0
                s03 = p3 - p0
                s21 = p1 - p2
                s23 = p3 - p2
                
                # First triangle
                t1norm = (s02.cross(s01)).norm()
                a = s01.mag()
                b = s02.mag()
                c = s21.mag()
                s = (a + b + c)/2
                A1 = np.sqrt(s*(s - a)*(s - b)*(s - c))/1000

                # Second triangle
                t2norm = (s03.cross(s02)).norm()
                a = s03.mag()
                b = s23.mag()
                c = s02.mag()
                s = (a + b + c)/2
                A2 = np.sqrt(s*(s - a)*(s - b)*(s - c))/1000

                # Up and North
                p1up = (p1up - p1).norm()
                p3up = (p3up - p3).norm()
                p1N = (p1N - p1).norm()
                p3N = (p3N - p3).norm()

                # Combine
                norm_vec = norm_vec + A1*t1norm + A2*t2norm
                north_vec = north_vec + A1*p1N + A2*p3N
                up_vec = up_vec + A1*p1up + A2*p3up

        norm_vec = norm_vec.norm()
        north_vec = north_vec.norm()
        up_vec = up_vec.norm()

        # Do I need to flip the vector because it is pointing down (i.e.,
        # right-hand rule is violated)?
        flip = np.sign(up_vec.dot(norm_vec))
        norm_vec = flip*norm_vec

        # Angle between up_vec and norm_vec is dip
        self._dip = np.arcsin(up_vec.cross(norm_vec).mag())*180/np.pi

        # Normal vector projected to horizontal plane
        nvph = (norm_vec - up_vec.dot(norm_vec)*up_vec).norm()

        # Dip direction is angle between nvph and north; strike is orthogonal.
        cp = nvph.cross(north_vec)
        sign = np.sign(cp.dot(up_vec))
        dp = nvph.dot(north_vec)
        strike = np.arctan2(sign*cp.mag(), dp)*180/np.pi - 90
        if strike < -180:
            strike = strike + 360
        self._strike = strike


    def getDepthToTop(self):
        """
        Returns:
            float: Depth to top of rupture in km.
        """
        return np.min(self._topdeps)


    def getQuadrilaterals(self):
        """
        Return a list of quadrilaterals. Unlike QuadRupture, these
        quadrilaterals are not restricted to be coplanar or have 
        horizontal top/bottom edges.
        
        Return:
            list: List of quadrilaterals, where each quadrilateral is
                a list of OQ Points.
        """
        ugroup = np.unique(self._group_index)
        ngroup = len(ugroup)
        qlist = []
        for i in range(ngroup):
            ind = np.where(self._group_index == ugroup[i])[0]
            nq = len(ind) - 1
            for j in range(nq):
                P0 = Point(self._toplons[j],
                           self._toplats[j],
                           self._topdeps[j])
                P1 = Point(self._toplons[j+1],
                           self._toplats[j+1],
                           self._topdeps[j+1])
                P2 = Point(self._botlons[j+1],
                           self._botlats[j+1],
                           self._botdeps[j+1])
                P3 = Point(self._botlons[j],
                           self._botlats[j],
                           self._botdeps[j])
                qlist.append([P0, P1, P2, P3])

        return qlist

    def computeRrup(self, lon, lat, depth):
        """
        Method for computing rupture distance. 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Rupture distance (km).

        """

        mesh_dx = self._mesh_dx

        #-----------------------------------------------------------------------
        # Sort out sites
        #-----------------------------------------------------------------------
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        x, y, z = latlon2ecef(lat, lon, depth)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))

        #-----------------------------------------------------------------------
        # Get mesh
        #-----------------------------------------------------------------------
        mx = []
        my = []
        mz = []
        u_groups = np.unique(self._group_index)
        n_groups = len(u_groups)
        for j in range(n_groups):
            g_ind = np.where(u_groups[j] == self._group_index)[0]
            nq = len(self._toplats[g_ind]) - 1
            for i in range(nq):
                q = [Point(self._toplons[g_ind[i]],
                           self._toplats[g_ind[i]],
                           self._topdeps[g_ind[i]]), 
                     Point(self._toplons[g_ind[i+1]],
                           self._toplats[g_ind[i+1]],
                           self._topdeps[g_ind[i+1]]),
                     Point(self._botlons[g_ind[i+1]],
                           self._botlats[g_ind[i+1]],
                           self._botdeps[g_ind[i+1]]),
                     Point(self._botlons[g_ind[i]],
                           self._botlats[g_ind[i]],
                           self._botdeps[g_ind[i]])
                     ]
                mesh = get_quad_mesh(q, dx = mesh_dx)
                mx.extend(list(np.reshape(mesh['x'], (-1,)) ))
                my.extend(list(np.reshape(mesh['y'], (-1,)) ))
                mz.extend(list(np.reshape(mesh['z'], (-1,)) ))
        mesh_mat = np.array([np.array(mx), np.array(my), np.array(mz)])

        #-----------------------------------------------------------------------
        # Compute distance
        #-----------------------------------------------------------------------
        dist = np.zeros_like(x)
        for i in range(len(x)):
            sitecol = sites_ecef[i,:].reshape([3, 1])
            dif = sitecol - mesh_mat
            distarray = np.sqrt(np.sum(dif * dif, axis=0))
            dist[i] = np.min(distarray)/1000.0 # convert to km

        dist = np.reshape(dist, oldshape)

        return dist


    def computeRjb(self, lon, lat, depth):
        """
        Method for computing Joyner-Boore distance. 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Joyner-Boore distance (km).

        """

        mesh_dx = self._mesh_dx

        #-----------------------------------------------------------------------
        # Sort out sites
        #-----------------------------------------------------------------------
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        x, y, z = latlon2ecef(lat, lon, depth)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))

        #-----------------------------------------------------------------------
        # Get mesh
        #-----------------------------------------------------------------------
        mx = []
        my = []
        mz = []
        u_groups = np.unique(self._group_index)
        n_groups = len(u_groups)
        for j in range(n_groups):
            g_ind = np.where(u_groups[j] == self._group_index)[0]
            nq = len(self._toplats[g_ind]) - 1
            for i in range(nq):
                q = [Point(self._toplons[g_ind[i]],
                           self._toplats[g_ind[i]],
                           0),
                     Point(self._toplons[g_ind[i+1]],
                           self._toplats[g_ind[i+1]],
                           0),
                     Point(self._botlons[g_ind[i+1]],
                           self._botlats[g_ind[i+1]],
                           0),
                     Point(self._botlons[g_ind[i]],
                           self._botlats[g_ind[i]],
                           0)
                     ]
                mesh = get_quad_mesh(q, dx = mesh_dx)
                mx.extend(list(np.reshape(mesh['x'], (-1,)) ))
                my.extend(list(np.reshape(mesh['y'], (-1,)) ))
                mz.extend(list(np.reshape(mesh['z'], (-1,)) ))
        mesh_mat = np.array([np.array(mx), np.array(my), np.array(mz)])

        #-----------------------------------------------------------------------
        # Compute distance
        #-----------------------------------------------------------------------
        dist = np.zeros_like(x)
        for i in range(len(x)):
            sitecol = sites_ecef[i,:].reshape([3, 1])
            dif = sitecol - mesh_mat
            distarray = np.sqrt(np.sum(dif * dif, axis=0))
            dist[i] = np.min(distarray)/1000.0 # convert to km

        dist = np.reshape(dist, oldshape)

        return dist

    def computeGC2(self, lon, lat, depth):
        """
        Method for computing version 2 of the Generalized Coordinate system 
        (GC2) by Spudich and Chiou OFR 2015-1028. 
        
        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            dict: Dictionary with keys for each of the GC2-related distances,
                which include 'rx', 'ry', 'ry0', 'U', 'T'. 
        """
        # This just hands off to the module-level method
        # NOTE: Not sure if the non-horizontal top edges of EdgeRupture will
        #       case problems. Should do some checking. It might be okay to
        #       bring quad vertices up to the surface in this case. 
        dict = _computeGC2(self, lon, lat, depth)
        return dict


class QuadRupture(Rupture):
    """
    Rupture class that represents the rupture surface as a combination of 
    quadrilaterals. Each quadrilateral must have horizontal top and bottom
    edges and must be coplanar. These restrictions make the computation of
    rupture distances more efficient. The number of points in the top edges
    must match the number of points in the bottom edge. 
    """

    def __init__(self, d, origin):
        """
        Create a QuadRupture instance from a GeoJSON dictionary and an Origin.

        Args:
           d (dict): Rupture GeoJSON dictionary.
           origin (Origin): Reference to a ShakeMap Origin object.

        Returns:
            QuadRupture instance.

        """

        polys = d['features'][0]['geometry']['coordinates'][0]
        n_polygons = len(polys)
        lon = []
        lat = []
        dep = []
        for i in range(n_polygons):
            p = polys[i]
            p_lons = [pt[0] for pt in p][0:-1]
            p_lats = [pt[1] for pt in p][0:-1]
            p_depths = [pt[2] for pt in p][0:-1]
            lon = lon + p_lons + [np.nan]
            lat = lat + p_lats + [np.nan]
            dep = dep + p_depths + [np.nan]

        self._geojson = d
        self._lon = lon
        self._lat = lat
        self._depth = dep
        self._origin = origin
        self._reference = d['features'][0]['properties']['reference']
        self._setQuadrilaterals()


    @classmethod
    def fromArrays(cls, lon, lat, depth, origin, reference = ''):
        """
        Class method for constructing a QuadRupture from arrays of longitude
        latitude, and depth. The arrays should start on the top and move in the 
        strike direction, then down to the bottom and move back (without 
        closing) followed by a nan. 


        Args:
            lon (array): Sequence of rupture longitude vertices in clockwise 
                order.
            lat (array): Sequence of rupture latitude vertices in clockwise 
                order.
            depth (array): Sequence of rupture depth vertices in clockwise order.
            origin (Origin): Reference to a ShakeMap Origin object.
            reference (str): String citeable reference for Rupture.

        """
        pass
        

#        return cls(d, origin)



    def getLength(self):
        """
        Compute length of rupture based on top edge in km.

        Returns:
            float: Length of rupture (km).

        """
        flength = 0
        for quad in self._quadrilaterals:
            flength = flength + get_quad_length(quad)
        return flength

    def getWidth(self):
        """
        Compute average rupture width (km) for all quadrilaterals defined for
        the rupture.

        Returns:
            float: Average width in km of all rupture quadrilaterals.
        """
        wsum = 0.0
        for quad in self._quadrilaterals:
            wsum = wsum + get_quad_width(quad)
        mwidth = (wsum / len(self._quadrilaterals)) / 1000.0
        return mwidth

    def getArea(self):
        """
        Compute area of rupture. 
        
        Returns: 
            float: Rupture area in square km. 

        """
        asum = 0.0
        for quad in self._quadrilaterals:
            w = get_quad_width(quad)
            l = get_quad_length(quad)
            asum = asum + w*l
        return asum
        

            

    @classmethod
    def fromTrace(cls, xp0, yp0, xp1, yp1, zp, widths, dips, origin,
                  strike=None, group_index = None, reference=""):
        """
        Create a QuadRupture instance from a set of vertices that define the top
        of the rupture, and an array of widths/dips.

        Each rupture quadrilaterial is defined by specifying the latitude,
        longitude, and depth of the two vertices on the top edges, which must
        have the dame depths. The other verticies are then constructed from
        the top edges and the width and dip of the quadrilateral. 

        Args:
            xp0 (array): Array or list of longitudes (floats) of p0.
            yp0 (array): Array or list of latitudes (floats) of p0.
            xp1 (array): Array or list of longitudes (floats) of p1.
            yp1 (array): Array or list of latitudes (floats) of p1.
            zp (array): Array or list of depths for each of the top of rupture
                rectangles (km).
            widths (array): Array of widths for each of rectangle (km).
            dips (array): Array of dips for each of rectangle (degrees).
            origin (Origin): Reference to a ShakeMap origin object.
            strike (array): If None then strike is computed from verticies of
                top edge of each quadrilateral. If a scalar, then all
                quadrilaterals are constructed assuming this strike direction.
                If an array with the same length as the trace coordinates then
                it specifies the strike for each quadrilateral.
            group_index (list): List of integers to indicate group index. If
                None then each quadrilateral is assumed to be in a different
                group since there is no guarantee that any of them are
                continuous.
            reference (str): String explaining where the rupture definition came
                from (publication style reference, etc.).

        Returns:
            QuadRupture instance.

        """
        if len(xp0) == len(yp0) == len(xp1) == len(
                yp1) == len(zp) == len(dips) == len(widths):
            pass
        else:
            raise ShakeMapException(
                'Number of xp0,yp0,xp1,yp1,zp,widths,dips points must be equal.')
        if strike is None:
            pass
        else:
            if (len(xp0) == len(strike)) | (len(strike) == 1):
                pass
            else:
                raise ShakeMapException(
                    'Strike must be None, scalar, or same length as '\
                    'trace coordinates.')

        if group_index is None:
            group_index = np.array(range(len(xp0)))

        # Convert dips to radians
        dips = np.radians(dips)

        # Ensure that all input sequences are numpy arrays
        xp0 = np.array(xp0, dtype='d')
        xp1 = np.array(xp1, dtype='d')
        yp0 = np.array(yp0, dtype='d')
        yp1 = np.array(yp1, dtype='d')
        zp = np.array(zp, dtype='d')
        widths = np.array(widths, dtype='d')
        dips = np.array(dips, dtype='d')

        # Get a projection object
        west = np.min((xp0.min(), xp1.min()))
        east = np.max((xp0.max(), xp1.max()))
        south = np.min((yp0.min(), yp1.min()))
        north = np.max((yp0.max(), yp1.max()))

        # Projected coordinates are in km
        proj = get_orthographic_projection(west, east, north, south)
        xp2 = np.zeros_like(xp0)
        xp3 = np.zeros_like(xp0)
        yp2 = np.zeros_like(xp0)
        yp3 = np.zeros_like(xp0)
        zpdown = np.zeros_like(zp)
        for i in range(0, len(xp0)):
            # Project the top edge coordinates
            p0x, p0y = proj(xp0[i], yp0[i])
            p1x, p1y = proj(xp1[i], yp1[i])

            # Get the rotation angle defined by these two points
            if strike is None:
                dx = p1x - p0x
                dy = p1y - p0y
                theta = np.arctan2(dx, dy)  # theta is angle from north
            elif len(strike) == 1:
                theta = np.radians(strike[0])
            else:
                theta = np.radians(strike[i])

            R = np.array([[np.cos(theta), -np.sin(theta)],
                          [np.sin(theta), np.cos(theta)]])

            # Rotate the top edge points into a new coordinate system (vertical
            # line)
            p0 = np.array([p0x, p0y])
            p1 = np.array([p1x, p1y])
            p0p = np.dot(R, p0)
            p1p = np.dot(R, p1)

            # Get right side coordinates in project, rotated system
            dz = np.sin(dips[i]) * widths[i]
            dx = np.cos(dips[i]) * widths[i]
            p3xp = p0p[0] + dx
            p3yp = p0p[1]
            p2xp = p1p[0] + dx
            p2yp = p1p[1]

            # Get right side coordinates in un-rotated projected system
            p3p = np.array([p3xp, p3yp])
            p2p = np.array([p2xp, p2yp])
            Rback = np.array([[np.cos(-theta), -np.sin(-theta)],
                              [np.sin(-theta), np.cos(-theta)]])
            p3 = np.dot(Rback, p3p)
            p2 = np.dot(Rback, p2p)
            p3x = np.array([p3[0]])
            p3y = np.array([p3[1]])
            p2x = np.array([p2[0]])
            p2y = np.array([p2[1]])

            # project lower edge points back to lat/lon coordinates
            lon3, lat3 = proj(p3x, p3y, reverse=True)
            lon2, lat2 = proj(p2x, p2y, reverse=True)

            xp2[i] = lon2
            xp3[i] = lon3
            yp2[i] = lat2
            yp3[i] = lat3
            zpdown[i] = zp[i] + dz

        #-----------------------------------------------------------------------
        # Create GeoJSON object
        #-----------------------------------------------------------------------

        coords = []
        u_groups = np.unique(group_index)
        n_groups = len(u_groups)
        for i in range(n_groups):
            ind = np.where(u_groups[i] == group_index)[0]
            lons = np.concatenate([ xp0[ind[0]].reshape((1,)), xp1[ind], 
                                    xp2[ind][::-1], xp3[ind][::-1][-1].reshape((1,)), 
                                    xp0[ind[0]].reshape((1,)) ])
            lats = np.concatenate([ yp0[ind[0]].reshape((1,)), yp1[ind], 
                                    yp2[ind][::-1], yp3[ind][::-1][-1].reshape((1,)), 
                                    yp0[ind[0]].reshape((1,)) ])
            deps = np.concatenate([ zp[ind[0]].reshape((1,)), zp[ind],
                                    zpdown[ind][::-1], zpdown[ind][::-1][-1].reshape((1,)),
                                    zp[ind[0]].reshape((1,)) ])
            
            poly = []
            for lon,lat,dep in zip(lons, lats, deps):
                poly.append([lon, lat, dep])
            coords.append(poly)
            
        d = {"type":"FeatureCollection", 
             "metadata":{},
             "features":[{
                    "type":"Feature",
                    "properties":{
                        "rupture type":"rupture extent",
                        "reference":reference,
                        },
                    "geometry":{
                        "type":"MultiPolygon",
                        "coordinates":[coords]
                        }
            }]}

        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, ShakeDateTime):
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                d['metadata'][k] = v

        return cls(d, origin)

    def writeTextFile(self, rupturefile):
        """
        Write rupture data to rupture file format as defined in ShakeMap 
        Software Guide.

        Note that this currently treats each quadrilateral as a separate
        polygon. This needs to be udpated. 

        Args:
            rupturefile (str): Filename of output data file OR file-like object.

        """
        if not hasattr(rupturefile, 'read'):
            f = open(rupturefile, 'wt')
        else:
            f = rupturefile  # just a reference to the input file-like object
        f.write('#%s\n' % self._reference)
        for quad in self.getQuadrilaterals():
            P0, P1, P2, P3 = quad
            f.write('%.4f %.4f %.4f\n' % (P0.latitude, P0.longitude, P0.depth))
            f.write('%.4f %.4f %.4f\n' % (P1.latitude, P1.longitude, P1.depth))
            f.write('%.4f %.4f %.4f\n' % (P2.latitude, P2.longitude, P2.depth))
            f.write('%.4f %.4f %.4f\n' % (P3.latitude, P3.longitude, P3.depth))
            f.write('%.4f %.4f %.4f\n' % (P0.latitude, P0.longitude, P0.depth))
            f.write(u'>\n')
        if not hasattr(rupturefile, 'read'):
            f.close()


    @classmethod
    def fromVertices(cls,
                     xp0, yp0, zp0, xp1, yp1, zp1,
                     xp2, yp2, zp2, xp3, yp3, zp3,
                     origin,
                     group_index = None,
                     reference=None):
        """
        Create a QuadDrupture instance from the vector of vertices that fully
        define the quadrilaterals. The points p0, ..., p3 are labeled below for
        a trapezoid:

        ::

              p0--------p1
             /          |
            /           |
           p3-----------p2

        All of the following vector arguments must have the same length.

        Args: 
            xp0 (array): Array or list of longitudes (floats) of p0.
            yp0 (array): Array or list of latitudes (floats) of p0.
            zp0 (array): Array or list of depths (floats) of p0.
            xp1 (array): Array or list of longitudes (floats) of p1.
            yp1 (array): Array or list of latitudes (floats) of p1.
            zp1 (array): Array or list of depths (floats) of p1.
            xp2 (array): Array or list of longitudes (floats) of p2.
            yp2 (array): Array or list of latitudes (floats) of p2.
            zp2 (array): Array or list of depths (floats) of p2.
            xp3 (array): Array or list of longitudes (floats) of p3.
            yp3 (array): Array or list of latitudes (floats) of p3.
            zp3 (array): Array or list of depths (floats) of p3.
            origin (Origin): Reference to a ShakeMap Origin object.
            group_index (list): List of integers to indicate group index. If
                None then each quadrilateral is assumed to be in a different
                group since there is no guarantee that any of them are 
                continuous.
            reference (str): String explaining where the rupture definition came
                from (publication style reference, etc.)

        Returns:
            QuadRupture object, where the rupture is modeled as a series of
                trapezoids.

        """
        if len(xp0) == len(yp0) == len(zp0) == len(xp1) == len(yp1) == len(zp1) == \
           len(xp2) == len(yp2) == len(zp2) == len(xp3) == len(yp3) == len(zp3):
            pass
        else:
            raise ShakeMapException('All vectors specifying quadrilateral '\
                                    'vertices must have the same length.')

        nq = len(xp0)
        if group_index is not None:
            if len(group_index) != nq:
                raise Exception("group_index must have same length as vertices.")
        else:
            group_index = np.array(range(nq))

        xp0 = np.array(xp0, dtype='d')
        yp0 = np.array(yp0, dtype='d')
        zp0 = np.array(zp0, dtype='d')
        xp1 = np.array(xp1, dtype='d')
        yp1 = np.array(yp1, dtype='d')
        zp1 = np.array(zp1, dtype='d')
        xp2 = np.array(xp2, dtype='d')
        yp2 = np.array(yp2, dtype='d')
        zp2 = np.array(zp2, dtype='d')
        xp3 = np.array(xp3, dtype='d')
        yp3 = np.array(yp3, dtype='d')
        zp3 = np.array(zp3, dtype='d')

        #-----------------------------------------------------------------------
        # Create GeoJSON object
        #-----------------------------------------------------------------------

        coords = []
        u_groups = np.unique(group_index)
        n_groups = len(u_groups)
        for i in range(n_groups):
            ind = np.where(u_groups[i] == group_index)[0]
            lons = np.concatenate([ xp0[ind[0]].reshape((1,)), xp1[ind],
                                    xp2[ind][::-1], xp3[ind][::-1][-1].reshape((1,)),
                                    xp0[ind[0]].reshape((1,)) ])
            lats = np.concatenate([ yp0[ind[0]].reshape((1,)), yp1[ind],
                                    yp2[ind][::-1], yp3[ind][::-1][-1].reshape((1,)),
                                    yp0[ind[0]].reshape((1,)) ])
            deps = np.concatenate([ zp0[ind[0]].reshape((1,)), zp1[ind],
                                    zp2[ind][::-1], zp3[ind][::-1][-1].reshape((1,)),
                                    zp0[ind[0]].reshape((1,)) ])

            poly = []
            for lon,lat,dep in zip(lons, lats, deps):
                poly.append([lon, lat, dep])
            coords.append(poly)


        d = {"type":"FeatureCollection", 
             "metadata":{},
             "features":[{
                    "type":"Feature",
                    "properties":{
                        "rupture type":"rupture extent",
                        "reference":reference,
                        },
                    "geometry":{
                        "type":"MultiPolygon",
                        "coordinates":[coords]
                        }
            }]}

        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, ShakeDateTime):
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                d['metadata'][k] = v
        if hasattr(origin, 'id'):
            d['metadata']['eventid'] = origin.id

        return cls(d, origin)



    def getQuadrilaterals(self):
        """
        Return a list of quadrilaterals.

        Returns:
            list: List of quadrilaterals where each quad is a tuple of four
                `Point <https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py>`__
                objects.
        """
        return copy.deepcopy(self._quadrilaterals)

    def getStrike(self):
        """
        Return strike angle. If rupture consists of multiple quadrilaterals, the
        average strike angle, weighted by quad length, is returned.
        Note: for ruptures with quads where the strike angle changes by 180 deg
        due to reverses in dip direction are problematic and not handeled well
        by this algorithm.

        Returns:
            float: Strike angle in degrees.

        """
        nq = len(self._quadrilaterals)
        strikes = np.zeros(nq)
        lengths = np.zeros(nq)
        for i in range(nq):
            P0 = self._quadrilaterals[i][0]
            P1 = self._quadrilaterals[i][1]
            strikes[i] = P0.azimuth(P1)
            lengths[i] = get_quad_length(self._quadrilaterals[i])
        x = np.sin(np.radians(strikes))
        y = np.cos(np.radians(strikes))
        xbar = np.sum(x * lengths) / np.sum(lengths)
        ybar = np.sum(y * lengths) / np.sum(lengths)
        return np.degrees(np.arctan2(xbar, ybar))

    def getDepthToTop(self):
        """
        Determine shallowest vertex of entire rupture.

        :returns:
            Shallowest depth of all vertices (float).
        """
        mindep = 9999999
        for quad in self._quadrilaterals:
            P0, P1, P2, P3 = quad
            depths = np.array([P0.depth, P1.depth, P2.depth, P3.depth])
            if np.min(depths) < mindep:
                mindep = np.min(depths)
        return mindep

    def getDip(self):
        """
        Return average dip of all quadrilaterals in the rupture.

        Returns:
           float: Average dip in degrees.

        """
        dipsum = 0.0
        for quad in self._quadrilaterals:
            N = get_quad_normal(quad)
            V = get_vertical_vector(quad)
            dipsum = dipsum + np.degrees(np.arccos(Vector.dot(N, V)))
        dip = dipsum / len(self._quadrilaterals)
        return dip


    def getIndividualWidths(self):
        """
        Return an array of rupture widths (km), one for each quadrilateral
        defined for the rupture.

        Returns:
            Array of quad widths in km of all rupture quadrilaterals.
        """
        nquad = self.getNumQuads()
        widths = np.zeros(nquad)
        for i in range(nquad):
            q = self._quadrilaterals[i]
            widths[i] = get_quad_width(q) / 1000.0
        return widths

    def getIndividualTopLengths(self):
        """
        Return an array of rupture lengths along top edge (km),
        one for each quadrilateral defined for the rupture.

        :returns:
            Array of lengths in km of top edge of quadrilaterals.
        """
        nquad = self.getNumQuads()
        lengths = np.zeros(nquad)
        for i in range(nquad):
            P0, P1, P2, P3 = self._quadrilaterals[i]
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            lengths[i] = (p1 - p0).mag() / 1000.0
        return lengths





    @staticmethod
    def _fixStrikeDirection(quad):
        P0, P1, P2, P3 = quad
        eps = 1e-6
        p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
        p1 = Vector.fromPoint(P1)
        p2 = Vector.fromPoint(P2)
        p1p0 = p1 - p0
        p2p0 = p2 - p0
        qnv = Vector.cross(p2p0, p1p0).norm()
        tmp = p0 + qnv
        tmplat, tmplon, tmpz = ecef2latlon(tmp.x, tmp.y, tmp.z)
        if (tmpz - P0.depth) < eps: # If True then do nothing
            fixed = quad
        else:
            newP0 = copy.deepcopy(P1)
            newP1 = copy.deepcopy(P0)
            newP2 = copy.deepcopy(P3)
            newP3 = copy.deepcopy(P2)
            fixed = [newP0, newP1, newP2, newP3]
        return fixed
        

    def _setQuadrilaterals(self):
        """
        Create internal list of N quadrilaterals. Reverses quad if dip direction is
        incorrect.
        """

        # Make sure arrays are numpy arrays. 
        self._lon = np.array(self._lon)
        self._lat = np.array(self._lat)
        self._depth = np.array(self._depth)

        # Find the nans, which tells is where the separate polygons/groups are
        group_ends = np.where(np.isnan(self._lon))[0]
        n_groups = len(group_ends)

        # Check that arrays are the same length
        if len(self._lon) != len(self._lat) != len(self._depth):
            raise IndexError(
                'Length of input lon, lat, depth arrays must be equal')

        # Construct quads
        group_start = 0

        self._quadrilaterals = []
        self._group_index = []
        groupind = 0
        for i in range(n_groups):
            lonseg = self._lon[group_start:group_ends[i]]
            latseg = self._lat[group_start:group_ends[i]]
            depthseg = self._depth[group_start:group_ends[i]]
            
            # Each group can have many contiguous quadrilaterals defined in it
            # separations (nans) between segments mean that segments are not
            # contiguous.

            npoints = len(lonseg)
            nquads = int((npoints - 4) / 2) + 1
            quad_start = 0
            quad_end = -1
            for j in range(nquads):
                P0 = Point(lonseg[quad_start],
                           latseg[quad_start],
                           depthseg[quad_start])
                P1 = Point(lonseg[quad_start + 1],
                           latseg[quad_start + 1],
                           depthseg[quad_start + 1])
                P2 = Point(lonseg[quad_end - 1],
                           latseg[quad_end - 1],
                           depthseg[quad_end - 1])
                P3 = Point(lonseg[quad_end],
                           latseg[quad_end],
                           depthseg[quad_end])
                quad = [P0, P1, P2, P3]

                # Enforce plane by moving P2 -- already close because of check
                # in read_rupture_file/is_quadrupture_class/is_quad
                
                dummy, fixed_quad = is_quad(quad)

                # Reverse quad if necessary
                fixed_quad = self._fixStrikeDirection(fixed_quad)

                self._quadrilaterals.append(fixed_quad)
                
                quad_start = quad_start + 1
                quad_end = quad_end - 1

            group_start = group_ends[i] + 1
            self._group_index.extend([groupind] * nquads)
            groupind = groupind + 1

    def _getGroupIndex(self):
        """
        Return a list of segment group indexes.

        Returns:
            list: Segment group indexes; length equals the number of quadrilaterals.
        """
        return copy.deepcopy(self._group_index)

    @property
    def lats(self):
        """
        Return an array of latitudes for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected 
        segments.

        Returns:
            array: Numpy array of closed-loop latitude values; disconnected
                segments are separated by nans. 
        """
        lats = []
        quads = self.getQuadrilaterals()
        groups = self._getGroupIndex()
        u_groups = np.unique(groups)
        ng = len(u_groups)
        for i in range(ng):
            q_ind = np.where(groups == u_groups[i])[0]
            nq = len(q_ind)
            top_lats = []
            bot_lats = []
            for j in range(nq):
                if j == 0:
                    top0 = [quads[q_ind[j]][0].latitude]
                    bot0 = [quads[q_ind[j]][3].latitude]
                    top_lats = top_lats + top0
                    bot_lats = bot_lats + bot0
                top_lats = top_lats + [quads[q_ind[j]][1].latitude]
                bot_lats = bot_lats + [quads[q_ind[j]][2].latitude]
            lats = lats + top_lats + bot_lats[::-1] + top0 + [np.nan]
        
        return np.array(lats)

    @property
    def lons(self):
        """
        Return an array of longitudes for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected 
        segments.

        Returns:
            array: Numpy array of closed-loop longitude values; disconnected
                segments are separated by nans. 
        """
        lons = []
        quads = self.getQuadrilaterals()
        groups = self._getGroupIndex()
        u_groups = np.unique(groups)
        ng = len(u_groups)
        for i in range(ng):
            q_ind = np.where(groups == u_groups[i])[0]
            nq = len(q_ind)
            top_lons = []
            bot_lons = []
            for j in range(nq):
                if j == 0:
                    top0 = [quads[q_ind[j]][0].longitude]
                    bot0 = [quads[q_ind[j]][3].longitude]
                    top_lons = top_lons + top0
                    bot_lons = bot_lons + bot0
                top_lons = top_lons + [quads[q_ind[j]][1].longitude]
                bot_lons = bot_lons + [quads[q_ind[j]][2].longitude]
            lons = lons + top_lons + bot_lons[::-1] + top0 + [np.nan]
        return np.array(lons)

    @property
    def depths(self):
        """
        Return an array of depths for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected 
        segments.

        Returns:
            array: Numpy array of closed-loop depths; disconnected
                segments are separated by nans. 
        """
        deps = []
        quads = self.getQuadrilaterals()
        groups = self._getGroupIndex()
        u_groups = np.unique(groups)
        ng = len(u_groups)
        for i in range(ng):
            q_ind = np.where(groups == u_groups[i])[0]
            nq = len(q_ind)
            top_deps = []
            bot_deps = []
            for j in range(nq):
                if j == 0:
                    top0 = [quads[q_ind[j]][0].depitude]
                    bot0 = [quads[q_ind[j]][3].depitude]
                    top_deps = top_deps + top0
                    bot_deps = bot_deps + bot0
                top_deps = top_deps + [quads[q_ind[j]][1].depitude]
                bot_deps = bot_deps + [quads[q_ind[j]][2].depitude]
            deps = deps + top_deps + bot_deps[::-1] + top0 + [np.nan]
        
        return np.array(deps)

    def getDeps(self):
        """
        Return a copy of the array of depths for the rupture verticies.

        Returns:
            array: Numpy array of latitude values.
        """
        return self._depth.copy()

    def getNumGroups(self):
        """
        Return a count of the number of rupture groups.

        Returns:
            int:Rnumber of rupture groups. 

        """
        return len(np.unique(self._group_index))

    def getNumQuads(self):
        """
        Return a count of the number of rupture quadrilaterals.

        Returns:
            int: Number of rupture quadrilaterals.
        """
        return len(self._quadrilaterals)

    def getRuptureAsArrays(self):
        """
        Return a 3-tuple of numpy arrays indicating X, Y, Z (lon,lat,depth)
        coordinates. Rupture groups are separated by numpy.NaN values.

        Returns:
            tuple: 3-tuple of numpy arrays indicating X,Y,Z (lon,lat,depth) 
                coordinates.
        """
        return (np.array(self._lon), np.array(self._lat), np.array(self._depth))

    def getRuptureAsMesh(self):
        """
        Return rupture segments as a OQ-Hazardlib Mesh object.

        Returns:
            Mesh (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
        """
        rupture = Mesh(self._lon, self._lat, self._depth)
        return rupture

    def computeRjb(self, lon, lat, depth):
        """
        Method for computing Joyner-Boore distance. 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Joyner-Boore distance (km).

        """

        #-----------------------------------------------------------------------
        # Sort out sites
        #-----------------------------------------------------------------------
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        x, y, z = latlon2ecef(lat, lon, depth)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))
        
        minrjb = np.ones(newshape, dtype=lon.dtype) * 1e16
        quads = self.getQuadrilaterals()

        for i in range(len(quads)):
            P0, P1, P2, P3 = quads[i]
            S0 = copy.deepcopy(P0)
            S1 = copy.deepcopy(P1)
            S2 = copy.deepcopy(P2)
            S3 = copy.deepcopy(P3)
            S0.depth = 0.0
            S1.depth = 0.0
            S2.depth = 0.0
            S3.depth = 0.0
            squad = [S0, S1, S2, S3]
            rjbdist = _quad_distance(squad, sites_ecef)
            minrjb = np.minimum(minrjb, rjbdist)

        minrjb = minrjb.reshape(oldshape)
        return minrjb


    def computeRrup(self, lon, lat, depth):
        """
        Method for computing rupture distance. 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Rupture distance (km).

        """

        #-----------------------------------------------------------------------
        # Sort out sites
        #-----------------------------------------------------------------------
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        x, y, z = latlon2ecef(lat, lon, depth)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))
        
        minrrup = np.ones(newshape, dtype=lon.dtype) * 1e16
        quads = self.getQuadrilaterals()

        for i in range(len(quads)):
            rrupdist = _quad_distance(quads[i], sites_ecef)
            minrrup = np.minimum(minrrup, rrupdist)

        minrrup = minrrup.reshape(oldshape)
        return minrrup

    def computeGC2(self, lon, lat, depth):
        """
        Method for computing version 2 of the Generalized Coordinate system 
        (GC2) by Spudich and Chiou OFR 2015-1028. 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            dict: Dictionary with keys for each of the GC2-related distances,
                which include 'rx', 'ry', 'ry0', 'U', 'T'. 
        """
        # This just hands off to the module-level method
        dict = _computeGC2(self, lon, lat, depth)
        return dict



def _computeGC2(rupture, lon, lat, depth):
    """
    Method for computing GC2 from a ShakeMap Rupture instance.

    Args:
        rupture (Rupture): ShakeMap rupture object.
        lon (array): Numpy array of site longitudes.
        lat (array): Numpy array of site latitudes.
        depth (array): Numpy array of site depths.

    Returns: 
        dict: Dictionary of GC2 distances. Keys include "T", "U", "rx" 
            "ry", "ry0". 
    """

    quadlist = rupture.getQuadrilaterals()
    quadgc2 = copy.deepcopy(quadlist)

    oldshape = lon.shape

    if len(oldshape) == 2:
        newshape = (oldshape[0] * oldshape[1], 1)
    else:
        newshape = (oldshape[0], 1)

    #---------------------------------------------------------------------------
    # Define a projection that spans sites and rupture
    #---------------------------------------------------------------------------

    all_lat = np.append(lat, rupture.lats)
    all_lon = np.append(lon, rupture.lons)

    west = np.nanmin(all_lon)
    east = np.nanmax(all_lon)
    south = np.nanmin(all_lat)
    north = np.nanmax(all_lat)
    proj = get_orthographic_projection(west, east, north, south)

    totweight = np.zeros(newshape, dtype=lon.dtype)
    GC2T = np.zeros(newshape, dtype=lon.dtype)
    GC2U = np.zeros(newshape, dtype=lon.dtype)

    #---------------------------------------------------------------------------
    # First sort out strike discordance and nominal strike prior to
    # starting the loop if there is more than one group/trace.
    #---------------------------------------------------------------------------
    group_ind = rupture._getGroupIndex()

    # Need group_ind as numpy array for sensible indexing...
    group_ind_np = np.array(group_ind)
    uind = np.unique(group_ind_np)
    n_groups = len(uind)

    if n_groups > 1:
        #-----------------------------------------------------------------------
        # The first thing we need to worry about is finding the coordinate
        # shift. U's origin is "selected from the two endpoints most
        # distant from each other." 
        #-----------------------------------------------------------------------

        # Need to get index of first and last quad
        # for each segment
        iq0 = np.zeros(n_groups, dtype='int16')
        iq1 = np.zeros(n_groups, dtype='int16')
        for k in uind:
            ii = [i for i, j in enumerate(group_ind) if j == uind[k]]
            iq0[k] = int(np.min(ii))
            iq1[k] = int(np.max(ii))

        #-----------------------------------------------------------------------
        # This is an iterator for each possible combination of traces
        # including trace orientations (i.e., flipped). 
        #-----------------------------------------------------------------------

        it_seg = it.product(it.combinations(uind, 2),
                            it.product([0, 1], [0, 1]))

        # Placeholder for the trace pair/orientation that gives the
        # largest distance. 
        dist_save = 0

        for k in it_seg:
            s0ind = k[0][0]
            s1ind = k[0][1]
            p0ind = k[1][0]
            p1ind = k[1][1]
            if p0ind == 0:
                P0 = quadlist[iq0[s0ind]][0]
            else:
                P0 = quadlist[iq1[s0ind]][1]
            if p1ind == 0:
                P1 = quadlist[iq1[s1ind]][0]
            else:
                P1 = quadlist[iq0[s1ind]][1]

            dist = geodetic.distance(P0.longitude, P0.latitude, 0.0,
                                     P1.longitude, P1.latitude, 0.0)
            if dist > dist_save:
                dist_save = dist
                A0 = P0
                A1 = P1

        #-----------------------------------------------------------------------
        # A0 and A1 are the furthest two segment endpoints, but we still
        # need to sort out which one is the "origin".
        #-----------------------------------------------------------------------

        # This goofy while-loop is to adjust the side of the rupture where the
        # origin is located
        dummy = -1
        while dummy < 0:
            A0.depth = 0
            A1.depth = 0
            p_origin = Vector.fromPoint(A0)
            a0 = Vector.fromPoint(A0)
            a1 = Vector.fromPoint(A1)
            ahat = (a1 - a0).norm()

            # Loop over traces
            e_j = np.zeros(n_groups)
            b_prime = [None] * n_groups
            for j in range(n_groups):
                P0 = quadlist[iq0[j]][0]
                P1 = quadlist[iq1[j]][1]
                P0.depth = 0
                P1.depth = 0
                p0 = Vector.fromPoint(P0)
                p1 = Vector.fromPoint(P1)
                b_prime[j] = p1 - p0
                e_j[j] = ahat.dot(b_prime[j])
            E = np.sum(e_j)

            # List of discordancy
            dc = [np.sign(a) * np.sign(E) for a in e_j]
            b = Vector(0, 0, 0)
            for j in range(n_groups):
                b.x = b.x + b_prime[j].x * dc[j]
                b.y = b.y + b_prime[j].y * dc[j]
                b.z = b.z + b_prime[j].z * dc[j]
            bhat = b.norm()
            dummy = bhat.dot(ahat)
            if dummy < 0:
                tmpA0 = copy.deepcopy(A0)
                tmpA1 = copy.deepcopy(A1)
                A0 = tmpA1
                A1 = tmpA0

        #-----------------------------------------------------------------------
        # To fix discordancy, need to flip quads and rearrange
        # the order of quadgc2
        #-----------------------------------------------------------------------

        # 1) flip quads
        for i in range(len(quadgc2)):
            if dc[group_ind[i]] < 0:
                quadgc2[i] = reverse_quad(quadgc2[i])

        # 2) rearrange quadlist order
        qind = np.arange(len(quadgc2))
        for i in range(n_groups):
            qsel = qind[group_ind_np == uind[i]]
            if dc[i] < 0:
                qrev = qsel[::-1]
                qind[group_ind_np == uind[i]] = qrev

        quadgc2old = copy.deepcopy(quadgc2)
        for i in range(len(qind)):
            quadgc2[i] = quadgc2old[qind[i]]

        # End of if-statement for adjusting group discordancy

    s_i = 0.0
    l_i = np.zeros(len(quadgc2))

    for i in range(len(quadgc2)):
        G0, G1, G2, G3 = quadgc2[i]

        # Compute u_i and t_i for this quad
        t_i = __calc_t_i(G0, G1, lat, lon, proj)
        u_i = __calc_u_i(G0, G1, lat, lon, proj)

        # Quad length (top edge)
        l_i[i] = get_quad_length(quadgc2[i])

        #-----------------------------------------------------------------------
        # Weight of segment, three cases
        #-----------------------------------------------------------------------

        # Case 3: t_i == 0 and 0 <= u_i <= l_i
        w_i = np.zeros_like(t_i)

        # Case 1:
        ix = t_i != 0
        w_i[ix] = (1.0 / t_i[ix]) * (np.arctan((l_i[i] -
            u_i[ix]) / t_i[ix]) - np.arctan(-u_i[ix] / t_i[ix]))

        # Case 2:
        ix = (t_i == 0) & ((u_i < 0) | (u_i > l_i[i]))
        w_i[ix] = 1 / (u_i[ix] - l_i[i]) - 1 / u_i[ix]

        totweight = totweight + w_i
        GC2T = GC2T + w_i * t_i

        if n_groups == 1:
            GC2U = GC2U + w_i * (u_i + s_i)
        else:
            if i == 0:
                qind = np.array(range(len(quadgc2)))
                l_kj = 0
                s_ij_1 = 0
            else:
                l_kj = l_i[(group_ind_np == group_ind_np[i]) & (qind < i)]
                s_ij_1 = np.sum(l_kj)

            # First endpoint in the current 'group' (or 'trace' in GC2 terms)
            p1 = Vector.fromPoint(quadgc2[iq0[group_ind[i]]][0])
            s_ij_2 = (p1 - p_origin).dot(np.sign(E)*ahat) / 1000.0

            # Above is GC2N, for GC2T use:
            #s_ij_2 = (p1 - p_origin).dot(bhat) / 1000.0

            s_ij = s_ij_1 + s_ij_2
            GC2U = GC2U + w_i * (u_i + s_ij)

        s_i = s_i + l_i[i]

    GC2T = GC2T / totweight
    GC2U = GC2U / totweight

    # Dictionary for holding the distances
    distdict = dict()
    
    distdict['T'] = copy.deepcopy(GC2T).reshape(oldshape)
    distdict['U'] = copy.deepcopy(GC2U).reshape(oldshape)

    # Take care of Rx
    Rx = copy.deepcopy(GC2T)  # preserve sign (no absolute value)
    Rx = Rx.reshape(oldshape)
    distdict['rx'] = Rx

    # Ry
    Ry = GC2U - s_i / 2.0
    Ry = Ry.reshape(oldshape)
    distdict['ry'] = Ry

    # Ry0
    Ry0 = np.zeros_like(GC2U)
    ix = GC2U < 0
    Ry0[ix] = np.abs(GC2U[ix])
    if n_groups > 1:
        s_i = s_ij + l_i[-1]
    ix = GC2U > s_i
    Ry0[ix] = GC2U[ix] - s_i
    Ry0 = Ry0.reshape(oldshape)
    distdict['ry0'] = Ry0

    return distdict





def get_quad_width(q):
    """
    Return width of an individual planar trapezoid, where the p0-p1 distance
    represents the long side.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        float: Width of planar trapezoid.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p3 = Vector.fromPoint(P3)
    AB = p0 - p1
    AC = p0 - p3
    t1 = (AB.cross(AC).cross(AB)).norm()
    width = t1.dot(AC)

    return width

def get_quad_mesh(q, dx):
    """
    Create a mesh from a quadrilateal. 

    Args:
        q (list): A quadrilateral; list of four points. 
        dx (float):  Target dx in km; used to get nx and ny of mesh, but mesh
            snaps to edges of rupture so actual dx/dy will not actually equal this
            value in general.

    Returns:
        dict: Mesh dictionary, which includes numpy arrays:

        - llx: lower left x coordinate in ECEF coords.
        - lly: lower left y coordinate in ECEF coords.
        - llz: lower left z coordinate in ECEF coords.
        - ulx: upper left x coordinate in ECEF coords.
        - etc.

    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)

    # Get nx based on length of top edge, minimum allowed is 2
    toplen_km = get_quad_length(q)
    nx = int(np.max([round(toplen_km / dx, 0) + 1, 2]))

    # Get array of points along top and bottom edges
    xfac = np.linspace(0, 1, nx)
    topp = [p0 + (p1 - p0) * a for a in xfac]
    botp = [p3 + (p2 - p3) * a for a in xfac]

    # Get ny based on mean length of vectors connecting top and bottom points
    ylen_km = np.ones(nx)
    for i in range(nx):
        ylen_km[i] = (topp[i] - botp[i]).mag() / 1000
    ny = int(np.max([round(np.mean(ylen_km) / dx, 0) + 1, 2]))
    yfac = np.linspace(0, 1, ny)

    # Build mesh: dict of ny by nx arrays (x, y, z):
    mesh = {'x': np.zeros([ny, nx]), 'y': np.zeros(
        [ny, nx]), 'z': np.zeros([ny, nx])}
    for i in range(nx):
        mpts = [topp[i] + (botp[i] - topp[i]) * a for a in yfac]
        mesh['x'][:, i] = [a.x for a in mpts]
        mesh['y'][:, i] = [a.y for a in mpts]
        mesh['z'][:, i] = [a.z for a in mpts]

    # Make arrays of pixel corners
    mesh['llx'] = mesh['x'][1:, 0:-1]
    mesh['lrx'] = mesh['x'][1:, 1:]
    mesh['ulx'] = mesh['x'][0:-1, 0:-1]
    mesh['urx'] = mesh['x'][0:-1, 1:]
    mesh['lly'] = mesh['y'][1:, 0:-1]
    mesh['lry'] = mesh['y'][1:, 1:]
    mesh['uly'] = mesh['y'][0:-1, 0:-1]
    mesh['ury'] = mesh['y'][0:-1, 1:]
    mesh['llz'] = mesh['z'][1:, 0:-1]
    mesh['lrz'] = mesh['z'][1:, 1:]
    mesh['ulz'] = mesh['z'][0:-1, 0:-1]
    mesh['urz'] = mesh['z'][0:-1, 1:]
    mesh['cpx'] = np.zeros_like(mesh['llx'])
    mesh['cpy'] = np.zeros_like(mesh['llx'])
    mesh['cpz'] = np.zeros_like(mesh['llx'])

    # i and j are indices over subruptures
    ni, nj = mesh['llx'].shape
    for i in range(0, ni):
        for j in range(0, nj):
            # Rupture corner points
            pp0 = Vector(
                mesh['ulx'][i, j], mesh['uly'][i, j], mesh['ulz'][i, j])
            pp1 = Vector(
                mesh['urx'][i, j], mesh['ury'][i, j], mesh['urz'][i, j])
            pp2 = Vector(
                mesh['lrx'][i, j], mesh['lry'][i, j], mesh['lrz'][i, j])
            pp3 = Vector(
                mesh['llx'][i, j], mesh['lly'][i, j], mesh['llz'][i, j])
            # Find center of quad
            mp0 = pp0 + (pp1 - pp0) * 0.5
            mp1 = pp3 + (pp2 - pp3) * 0.5
            cp = mp0 + (mp1 - mp0) * 0.5
            mesh['cpx'][i, j] = cp.x
            mesh['cpy'][i, j] = cp.y
            mesh['cpz'][i, j] = cp.z
    return mesh


def get_local_unit_slip_vector(strike, dip, rake):
    """
    Compute the components of a unit slip vector.

    Args:
        strike (float): Clockwise angle (deg) from north of the line at the
            intersection of the rupture plane and the horizontal plane.
        dip (float): Angle (deg) between rupture plane and the horizontal plane
            normal to the strike (0-90 using right hand rule).
        rake (float): Direction of motion of the hanging wall relative to the
            foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in 'local' N-S, E-W, U-D coordinates.

    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.sin(rake) * np.cos(dip) * np.cos(strike) + \
        np.cos(rake) * np.sin(strike)
    sy = np.sin(rake) * np.cos(dip) * np.sin(strike) + \
        np.cos(rake) * np.cos(strike)
    sz = np.sin(rake) * np.sin(dip)
    return Vector(sx, sy, sz)


def get_local_unit_slip_vector_DS(strike, dip, rake):
    """
    Compute the DIP SLIP components of a unit slip vector.

    Args:
        strike (float): Clockwise angle (deg) from north of the line at the
            intersection of the rupture plane and the horizontal plane.
        dip (float): Angle (degrees) between rupture plane and the horizontal
            plane normal to the strike (0-90 using right hand rule).
        rake (float): Direction of motion of the hanging wall relative to the
            foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in 'local' N-S, E-W, U-D coordinates.

    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.sin(rake) * np.cos(dip) * np.cos(strike)
    sy = np.sin(rake) * np.cos(dip) * np.sin(strike)
    sz = np.sin(rake) * np.sin(dip)
    return Vector(sx, sy, sz)


def get_local_unit_slip_vector_SS(strike, dip, rake):
    """
    Compute the STRIKE SLIP components of a unit slip vector.

    Args:
        strike (float): Clockwise angle (deg) from north of the line at the
            intersection of the rupture plane and the horizontal plane.
        dip (float): Angle (degrees) between rupture plane and the horizontal
            plane normal to the strike (0-90 using right hand rule).
        rake (float): Direction of motion of the hanging wall relative to the
            foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in 'local' N-S, E-W, U-D coordinates.

    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.cos(rake) * np.sin(strike)
    sy = np.cos(rake) * np.cos(strike)
    sz = 0.0
    return Vector(sx, sy, sz)

def reverse_quad(q):
    """
    Reverse the verticies of a quad in the sense that the strike direction
    is flipped. 

    Args:
        q (list): A quadrilateral; list of four points.

    Returns: 
        list: Reversed quadrilateral.

    """
    return [q[1], q[0], q[3], q[2]]

def get_quad_slip(q, rake):
    """
    Compute the unit slip vector in ECEF space for a quad and rake angle.

    Args:
        q (list): A quadrilateral; list of four points.
        rake (float): Direction of motion of the hanging wall relative to
        the foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in ECEF space.

    """
    P0, P1, P2 = q[0:3]
    strike = P0.azimuth(P1)
    dip = get_quad_dip(q)
    s1_local = get_local_unit_slip_vector(strike, dip, rake)
    s0_local = Vector(0, 0, 0)
    qlats = [a.latitude for a in q]
    qlons = [a.longitude for a in q]
    proj = get_orthographic_projection(
        np.min(qlons), np.max(qlons), np.min(qlats), np.max(qlats))
    s1_ll = proj(np.array([s1_local.x]), np.array([s1_local.y]), reverse=True)
    s0_ll = proj(np.array([s0_local.x]), np.array([s0_local.y]), reverse=True)
    s1_ecef = Vector.fromTuple(latlon2ecef(s1_ll[1], s1_ll[0], s1_local.z))
    s0_ecef = Vector.fromTuple(latlon2ecef(s0_ll[1], s0_ll[0], s0_local.z))
    slp_ecef = (s1_ecef - s0_ecef).norm()
    return slp_ecef


def get_quad_length(q):
    """
    Length of top eduge of a quadrilateral.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        float: Length of quadrilateral in km.

    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    qlength = (p1 - p0).mag() / 1000
    return qlength


def get_quad_dip(q):
    """
    Dip of a quadrilateral.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        float: Dip in degrees.

    """
    N = get_quad_normal(q)
    V = get_vertical_vector(q)
    dip = np.degrees(np.arccos(Vector.dot(N, V)))
    return dip


def get_quad_normal(q):
    """
    Compute the unit normal vector for a quadrilateral in
    ECEF coordinates.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: Normalized normal vector for the quadrilateral in ECEF coords.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    p3 = Vector.fromPoint(P3)
    v1 = p1 - p0
    v2 = p3 - p0
    vn = Vector.cross(v2, v1).norm()
    return vn


def get_quad_strike_vector(q):
    """
    Compute the unit vector pointing in the direction of strike for a
    quadrilateral in ECEF coordinates. Top edge assumed to be horizontal.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: The unit vector pointing in strike direction in ECEF coords.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    v1 = (p1 - p0).norm()
    return v1


def get_quad_down_dip_vector(q):
    """
    Compute the unit vector pointing down dip for a quadrilateral in
    ECEF coordinates.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: The unit vector pointing down dip in ECEF coords.

    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    p0p1 = p1 - p0
    qnv = get_quad_normal(q)
    ddv = Vector.cross(p0p1, qnv).norm()
    return ddv


def get_vertical_vector(q):
    """
    Compute the vertical unit vector for a quadrilateral
    in ECEF coordinates.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: Normalized vertical vector for the quadrilateral in ECEF coords.
    """
    P0, P1, P2, P3 = q
    P0_up = copy.deepcopy(P0)
    P0_up.depth = P0_up.depth - 1.0
    p0 = Vector.fromPoint(P0)   # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P0_up)
    v1 = (p1 - p0).norm()
    return v1

def _quad_distance(q, points):
    """
    Calculate the shortest distance from a set of points to a rupture surface.

    Args:
        q (list): A quadrilateral; list of four points.
        points (array): Numpy array Nx3 of points (ECEF) to calculate distance 
            from.

    Returns:
        float: Array of size N of distances (in km) from input points to rupture
            surface.
    """
    P0, P1, P2, P3 = q

    # Convert to ecef
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)

    # Make a unit vector normal to the plane
    normalVector = (p1 - p0).cross(p2 - p0).norm()

    dist = np.ones_like(points[:, 0]) * np.nan

    p0d = p0.getArray() - points
    p1d = p1.getArray() - points
    p2d = p2.getArray() - points
    p3d = p3.getArray() - points

    # Create 4 planes with normals pointing outside rectangle
    n0 = (p1 - p0).cross(normalVector).getArray()
    n1 = (p2 - p1).cross(normalVector).getArray()
    n2 = (p3 - p2).cross(normalVector).getArray()
    n3 = (p0 - p3).cross(normalVector).getArray()

    sgn0 = np.signbit(np.sum(n0 * p0d, axis=1))
    sgn1 = np.signbit(np.sum(n1 * p1d, axis=1))
    sgn2 = np.signbit(np.sum(n2 * p2d, axis=1))
    sgn3 = np.signbit(np.sum(n3 * p3d, axis=1))

    inside_idx = (sgn0 == sgn1) & (sgn1 == sgn2) & (sgn2 == sgn3)
    dist[inside_idx] = np.power(np.abs(
        np.sum(p0d[inside_idx, :] * normalVector.getArray(), axis=1)), 2)

    outside_idx = np.logical_not(inside_idx)
    s0 = _distance_sq_to_segment(p0d, p1d)
    s1 = _distance_sq_to_segment(p1d, p2d)
    s2 = _distance_sq_to_segment(p2d, p3d)
    s3 = _distance_sq_to_segment(p3d, p0d)

    smin = np.minimum(np.minimum(s0, s1), np.minimum(s2, s3))
    dist[outside_idx] = smin[outside_idx]
    dist = np.sqrt(dist) / 1000.0
    shp = dist.shape
    if len(shp) == 1:
        dist.shape = (shp[0], 1)
    if np.any(np.isnan(dist)):
        raise ShakeMapException("Could not calculate some distances!")
    dist = np.fliplr(dist)
    return dist


def get_distance_to_plane(planepoints, otherpoint):
    """
    Calculate a point's distance to a plane.  Used to figure out if a
    quadrilateral points are all co-planar.

    Args:
        planepoints (list): List of three points (from Vector class) defining a
            plane.
        otherpoint (Vector): 4th Vector to compare to points defining the plane.
        
    Returns:
        float: Distance (in meters) from otherpoint to plane.

    """
    # from
    # https://en.wikipedia.org/wiki/Plane_(geometry)#Describing_a_plane_through_three_points
    p0, p1, p2 = planepoints
    x1, y1, z1 = p0.getArray()
    x2, y2, z2 = p1.getArray()
    x3, y3, z3 = p2.getArray()
    D = np.linalg.det(np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]))
    if D != 0:
        d = -1
        at = np.linalg.det(np.array([[1, y1, z1], [1, y2, z2], [1, y3, z3]]))
        bt = np.linalg.det(np.array([[x1, 1, z1], [x2, 1, z2], [x3, 1, z3]]))
        ct = np.linalg.det(np.array([[x1, y1, 1], [x2, y2, 1], [x3, y3, 1]]))
        a = (-d / D) * at
        b = (-d / D) * bt
        c = (-d / D) * ct

        numer = np.abs(a * otherpoint.x +
                       b * otherpoint.y +
                       c * otherpoint.z + d)
        denom = np.sqrt(a**2 + b**2 + c**2)
        dist = numer / denom
    else:
        dist = 0
    return dist

def __calc_u_i(P0, P1, lat, lon, proj):
    """
    Calculate u_i distance. See Spudich and Chiou OFR 2015-1028. This is the 
    distance along strike from the first vertex (P0) of the i-th segment.

    Args:
        P0 (point): OQ Point object, representing the first top-edge vertex of 
            a rupture quadrilateral.
        P1 (point): OQ Point object, representing the second top-edge vertex of 
            a rupture quadrilateral.
        lat (array): A numpy array of latitudes.
        lon (array): A numpy array of longitudes.
        proj (object): An orthographic projection from 
            openquake.hazardlib.geo.utils.get_orthographic_projection.

    Returns:
        array: Array of size N of distances (in km) from input points to rupture
            surface.

    """

    # projected coordinates are in km
    p0x, p0y = proj(P0.x, P0.y)
    p1x, p1y = proj(P1.x, P1.y)

    # Unit vector pointing along strike
    u_i_hat = Vector(p1x - p0x, p1y - p0y, 0).norm()

    # Convert sites to Cartesian
    sx, sy = proj(lon, lat)
    sx1d = np.reshape(sx, (-1,))
    sy1d = np.reshape(sy, (-1,))

    # Vectors from P0 to sites
    r = np.zeros([len(sx1d), 2])
    r[:, 0] = sx1d - p0x
    r[:, 1] = sy1d - p0y

    # Dot product gives u_i
    u_i = np.sum(u_i_hat.getArray()[0:2] * r, axis=1)
    shp = u_i.shape
    if len(shp) == 1:
        u_i.shape = (shp[0], 1)
    u_i = np.fliplr(u_i)

    return u_i

def __calc_t_i(P0, P1, lat, lon, proj):
    """
    Calculate t_i distance. See Spudich and Chiou OFR 2015-1028. This is the
    distance measured normal to strike from the i-th segment. Values on the
    hanging-wall are positive and those on the foot-wall are negative.

    Args:
        P0 (point): OQ Point object, representing the first top-edge vertex of 
            a rupture quadrilateral.
        P1 (point): OQ Point object, representing the second top-edge vertex of 
            a rupture quadrilateral.
        lat (array): A numpy array of latitudes.
        lon (array): A numpy array of longitudes.
        proj (object): An orthographic projection from 
            openquake.hazardlib.geo.utils.get_orthographic_projection.

    Returns:
        array: Array of size N of distances (in km) from input points to rupture
            surface.
    """

    # projected coordinates are in km
    p0x, p0y = proj(P0.x, P0.y)
    p1x, p1y = proj(P1.x, P1.y)

    # Unit vector pointing normal to strike
    t_i_hat = Vector(p1y - p0y, -(p1x - p0x), 0).norm()

    # Convert sites to Cartesian
    sx, sy = proj(lon, lat)
    sx1d = np.reshape(sx, (-1,))
    sy1d = np.reshape(sy, (-1,))

    # Vectors from P0 to sites
    r = np.zeros([len(sx1d), 2])
    r[:, 0] = sx1d - p0x
    r[:, 1] = sy1d - p0y

    # Dot product gives t_i
    t_i = np.sum(t_i_hat.getArray()[0:2] * r, axis=1)
    shp = t_i.shape
    if len(shp) == 1:
        t_i.shape = (shp[0], 1)
    t_i = np.fliplr(t_i)
    return t_i


def _distance_sq_to_segment(p0, p1):
    """
    Calculate the distance^2 from the origin to a segment defined by two vectors

    Args:
        p0 (array): Numpy array Nx3 (ECEF).
        p1 (array): Numpy array Nx3 (ECEF).
    
    Returns:
        float: The squared distance from the origin to a segment.
    """
    # /*
    #  * This algorithm is from (Vince's) CS1 class.
    #  * It returns the distance^2 from the origin to a segment defined
    #  * by two vectors
    #  */

    dist = np.zeros_like(p1[:, 0])
    # /* Are the two points equal? */
    idx_equal = (p0[:, 0] == p1[:, 0]) & (
        p0[:, 1] == p1[:, 1]) & (p0[:, 2] == p1[:, 2])
    dist[idx_equal] = np.sqrt(p0[idx_equal, 0]**2 +
                              p0[idx_equal, 1]**2 + p0[idx_equal, 2]**2)

    v = p1 - p0

    # /*
    #  * C1 = c1/|v| is the projection of the origin O on line (P0,P1).
    #  * If C1 is negative, then O is outside the segment and
    #  * closer to the P0 side.
    #  * If C1 is positive and >V then O is on the other side.
    #  * If C1 is positive and <V then O is inside.
    #  */

    c1 = -1 * np.sum(p0 * v, axis=1)
    idx_neg = c1 <= 0
    dist[idx_neg] = p0[idx_neg, 0]**2 + p0[idx_neg, 1]**2 + p0[idx_neg, 2]**2

    c2 = np.sum(v * v, axis=1)
    idx_less_c1 = c2 <= c1
    dist[idx_less_c1] = p1[idx_less_c1, 0]**2 + \
        p1[idx_less_c1, 1]**2 + p1[idx_less_c1, 2]**2

    idx_other = np.logical_not(idx_neg | idx_equal | idx_less_c1)

    nr, nc = p0.shape
    t1 = c1 / c2
    t1.shape = (nr, 1)
    tmp = p0 + (v * t1)
    dist[idx_other] = tmp[idx_other, 0]**2 + \
        tmp[idx_other, 1]**2 + tmp[idx_other, 2]**2

    return dist

def rake_to_mech(rake):
    """
    Convert rake to mechanism.

    Args:
        rake (float): Rake angle in degrees.

    Returns:
        str: Mechanism.
    """
    mech = 'ALL'
    if rake is not None:
        if (rake >= -180 and rake <= -150) or \
           (rake >= -30  and rake <= 30) or \
           (rake >= 150 and rake <= 180):
            mech = 'SS'
        if rake >= -120 and rake <= -60:
            mech = 'NM'
        if rake >= 60 and rake <= 120:
            mech = 'RS'
    return mech
