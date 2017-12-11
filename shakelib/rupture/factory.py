#!/usr/bin/env python3

# stdlib modules
import json

# third party imports
import numpy as np
from openquake.hazardlib.geo.point import Point

# local imports
from shakelib.rupture.point_rupture import PointRupture
from shakelib.rupture.quad_rupture import QuadRupture
from shakelib.rupture.edge_rupture import EdgeRupture
from shakelib.rupture.origin import Origin
from shakelib.rupture import utils
from shakelib.rupture import constants
from shakelib.utils.exception import ShakeLibException


def get_rupture(origin, file=None, mesh_dx=0.5):
    """
    This is a module-level function to read in a rupture file. This allows for
    the ShakeMap 3 text file specification or the ShakeMap 4 JSON rupture
    format. The ShakeMap 3 (".txt" extension) only supports QuadRupture style
    rupture representation and so this method will always return a QuadRupture
    instance. The ShakeMap 4 JSON format supports QuadRupture and EdgeRupture
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
            # -----------------------------------------------------------------
            # First, try to read as a json file
            # -----------------------------------------------------------------
            if isinstance(file, str):
                with open(file) as f:
                    d = json.load(f)
            else:
                d = json.loads(str(file))

            rupt = rupture_from_dict_and_origin(d, origin, mesh_dx=mesh_dx)

        except json.JSONDecodeError:
            # -----------------------------------------------------------------
            # Reading as json failed, so hopefully it is a ShakeMap 3 text file
            # -----------------------------------------------------------------
            try:
                d = text_to_json(file)
                rupt = rupture_from_dict_and_origin(d, origin, mesh_dx=mesh_dx)
            except:
                raise Exception("Unknown rupture file format.")
    else:
        if origin is None:
            raise Exception("Origin requred if no rupture file is provided.")
        rupt = PointRupture(origin)
    return rupt


def rupture_from_dict_and_origin(d, origin, mesh_dx=0.5):
    """
    Method returns either a QuadRupture or EdgeRupture object based on a
    GeoJSON dictionary and an origin. Note that this is very similar to
    :func:`rupture_from_dict`, except that method is for
    constructing the rupture objects from a dict that already contains the
    origin info in the `metadata` field (e.g., from a dict from a Shakemap
    container), while this method is for construction of the rupture objects
    from a GeoJSON dict that does not yet include that information (e.g., from
    a dict that is read in to initially create the shakemap container, along
    with an Origin that is derived from `event.xml`).

    .. seealso:: :func:`rupture_from_dict`

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
        rupt = EdgeRupture(d, origin, mesh_dx=mesh_dx)

    return rupt


def rupture_from_dict(d):
    """
    Method returns either a Rupture subclass (QuadRupture, EdgeRupture, or
    PointRupture) object based on a GeoJSON dictionary.

    .. seealso::
        :func:`rupture_from_dict_and_origin`

    Args:
        d (dict): Rupture GeoJSON dictionary, which must contain origin
            information in the 'metadata' field.

    Returns:
        a Rupture subclass.

    """
    validate_json(d)

    origin = Origin(d['metadata'])

    # What type of rupture is this?
    geo_type = d['features'][0]['geometry']['type']
    if geo_type == 'MultiPolygon':
        # EdgeRupture will have 'mesh_dx' in metadata
        if 'mesh_dx' in d['metadata']:
            mesh_dx = d['metadata']['mesh_dx']
            rupt = EdgeRupture(d, origin, mesh_dx=mesh_dx)
        else:
            rupt = QuadRupture(d, origin)
    elif geo_type == 'Point':
        rupt = PointRupture(origin)

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
                * Verticies within a rupture group must start along the top
                  edge and move in the strike direction then move to the bottom
                  edge and move back in the opposite direction.

    Returns:
        dict: GeoJSON rupture dictionary.

    """

    # -------------------------------------------------------------------------
    # First read in the data
    # -------------------------------------------------------------------------
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
            raise ShakeLibException(
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
        "metadata": {
            'reference':reference
        },
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "rupture type": "rupture extent"
                },
                "geometry": {
                    "type": "MultiPolygon",
                    "coordinates": [coords]
                }
            }
        ]
    }
    return d


def validate_json(d):
    """
    Check that the JSON format is acceptable. This is only for requirements
    that are common to both QuadRupture and EdgeRupture.

    Args:
        d (dict): Rupture JSON dictionary.
    """
    if d['type'] != 'FeatureCollection':
        raise Exception('JSON file is not a \"FeatureColleciton\".')

    if len(d['features']) != 1:
        raise Exception('JSON file should contain excactly one feature.')

    if 'reference' not in d['metadata'].keys():
        raise Exception('Json metadata field should contain '
                        '\"reference\" key.')

    f = d['features'][0]

    if f['type'] != 'Feature':
        raise Exception('Feature type should be \"Feature\".')

    geom = f['geometry']

    if (geom['type'] != 'MultiPolygon' and
            geom['type'] != 'Point'):
        raise Exception('Geometry type should be \"MultiPolygon\" '
                        'or \"Point\".')

    if 'coordinates' not in geom.keys():
        raise Exception('Geometry dictionary should contain \"coordinates\" '
                        'key.')

    polygons = geom['coordinates'][0]

    if geom['type'] == 'MultiPolygon':
        n_polygons = len(polygons)
        for i in range(n_polygons):
            p = polygons[i]
            n_points = len(p)
            if n_points % 2 == 0:
                raise Exception('Number of points in polyon must be odd.')

            if p[0] != p[-1]:
                raise Exception('First and last points in polygon must be '
                                'identical.')

            n_pairs = int((n_points - 1) / 2)
            for j in range(n_pairs):
                # -------------------------------------------------------------
                # Points are paired and in each pair the top is first, as in:
                #
                #      _.-P1-._
                #   P0'        'P2---P3
                #   |                  \
                #   P7---P6----P5-------P4
                #
                # Pairs: P0-P7, P1-P6, P2-P5, P3-P4
                # -------------------------------------------------------------
                top_depth = p[j][2]
                bot_depth = p[-(j + 2)][2]
                if top_depth >= bot_depth:
                    raise Exception(
                        'Top points must be ordered before bottom points.')


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
        n_pairs = int((n_points - 1) / 2)
        n_quads = n_pairs - 1

        for k in range(n_quads):
            # Four points of each quad should be co-planar within a tolerance
            quad = [Point(p[k][0], p[k][1], p[k][2]),
                    Point(p[k + 1][0], p[k + 1][1], p[k + 1][2]),
                    Point(p[-(k + 3)][0], p[-(k + 3)][1], p[-(k + 3)][2]),
                    Point(p[-(k + 2)][0], p[-(k + 2)][1], p[-(k + 2)][2])]
            test = utils.is_quad(quad)
            if test[0] is False:
                isQuad = False

            # Within each quad, top and bottom edges must be horizontal
            tops = np.array([quad[0].depth, quad[1].depth])
            if not np.isclose(tops[0], tops, rtol=0,
                              atol=constants.DEPTH_TOL).all():
                isQuad = False
            bots = np.array([quad[2].depth, quad[3].depth])
            if not np.isclose(bots[0], bots, rtol=0,
                              atol=constants.DEPTH_TOL).all():
                isQuad = False

    return isQuad
