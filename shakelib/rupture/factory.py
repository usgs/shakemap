#!/usr/bin/env python3

# stdlib modules
import json
import logging
import copy

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
from impactutils.time.ancient_time import HistoricTime


def get_rupture(origin, file=None, mesh_dx=0.5, new_format=True):
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
        origin (Origin):
            A ShakeMap origin instance; required because hypocentral/epicentral
            distances are computed from the Ruptureclass.

        file (srt):
            Path to rupture file (optional).

        mesh_dx (float):
            Target spacing (in km) for rupture discretization; default is
            0.5 km and it is only used if the rupture file is an EdgeRupture.

        new_format (bool):
            Indicates whether text rupture format is
            "old" (lat, lon, depth) or "new" (lon, lat, depth) style.

    Returns:
        Rupture subclass instance.

    """
    if not isinstance(origin, Origin):
        raise TypeError("An Origin is requred.")

    if file is not None:
        try:
            # -----------------------------------------------------------------
            # First, try to read as a json file
            # -----------------------------------------------------------------
            if isinstance(file, str):
                with open(file, encoding="latin-1") as f:
                    d = json.load(f)
            else:
                d = json.loads(str(file))

            rupt = rupture_from_dict_and_origin(d, origin, mesh_dx=mesh_dx)

        except Exception as e:
            if not isinstance(e, json.JSONDecodeError) and \
               not isinstance(e, UnicodeDecodeError):
                logging.warning("Unknown exception reading fault file: %s" %
                                str(e))
            # -----------------------------------------------------------------
            # Reading as json failed, so hopefully it is a ShakeMap 3 text file
            # -----------------------------------------------------------------
            try:
                d = text_to_json(file, new_format=new_format)
                rupt = rupture_from_dict_and_origin(d, origin, mesh_dx=mesh_dx)
            except ValueError as e:
                logging.error(e)
                raise ValueError("Error reading %s. This could have been "
                                 "because the latitude and longitdue are "
                                 "reversed. Try changing the 'new_format' "
                                 "option." % file)
            except Exception as e:
                logging.error(e)
                raise IOError("Unknown rupture file format.")
    else:
        rupt = PointRupture(origin)

    return rupt


def rupture_from_dict_and_origin(rupdict, origin, mesh_dx=0.5):
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
        rupdictd (dict):
            Rupture GeoJSON dictionary.
        origin (Origin):
            A ShakeMap origin object.
        mesh_dx (float):
            Target spacing (in km) for rupture discretization;
            default is 0.5 km and it is only used if the rupture file is an
            EdgeRupture.

    Returns:
        a Rupture subclass.

    """
    validate_json(rupdict)

    # Is this a QuadRupture or an EdgeRupture?
    valid_quads = is_quadrupture_class(rupdict)

    if valid_quads is True:
        rupt = QuadRupture(rupdict, origin)
    else:
        if rupdict['features'][0]['geometry']['type'] == 'Point':
            rupt = PointRupture(origin)
        else:
            rupt = EdgeRupture(rupdict, origin, mesh_dx=mesh_dx)

    return rupt


def rupture_from_dict(d):
    """
    Method returns either a Rupture subclass (QuadRupture, EdgeRupture, or
    PointRupture) object based on a GeoJSON dictionary.

    .. seealso::
        :func:`rupture_from_dict_and_origin`

    Args:
        d (dict):
            Rupture GeoJSON dictionary, which must contain origin
            information in the 'metadata' field.

    Returns:
        a Rupture subclass.

    """
    validate_json(d)

    # We don't want to mess with the input just in case it gets used again
    d = copy.deepcopy(d)

    try:
        d['metadata']['time'] = HistoricTime.strptime(d['metadata']['time'],
                                                      constants.TIMEFMT)
    except ValueError:
        d['metadata']['time'] = HistoricTime.strptime(d['metadata']['time'],
                                                      constants.ALT_TIMEFMT)

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


def _rotate_polygon(p):
    """
    This is a method to rotate a polygon, which is used to try to
    fix incorrectly specified polygons.

    Args:
        p (list):
            A list of lon/lat/depth lists.

    Returns:
        list: Same as before but the points have been incremented
        by one in position.

    """
    # Convert to numpy array for better indexing support
    parray = np.array(p)

    # Ensure the polygon is closed
    if not np.array_equal(parray[-1, :], parray[0, :]):
        raise ShakeLibException('Rupture file has unclosed segments.')

    # Drop last point
    parray = parray[0:-1, :]

    # Rotate by putting first at end
    parray = np.append(parray[1:, :],
                       np.reshape(parray[0, :], (1, -1)),
                       axis=0)

    # Put new first point onto the end so that it is closed.
    parray = np.append(parray,
                       np.reshape(parray[0, :], (1, -1)),
                       axis=0)

    # Turn numpy array back into a list
    polygon = parray.tolist()
    return polygon


def _check_polygon(p):
    """
    Check if the verticies are specified top first.

    Args:
        p (list):
            A list of five lon/lat/depth lists.

    Raises:
        ValueError: incorrectly specified polygon.

    """
    n_points = len(p)
    if n_points % 2 == 0:
        raise ValueError('Number of points in polyon must be odd.')

    if p[0] != p[-1]:
        raise ValueError('First and last points in polygon must be '
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
            raise ValueError(
                'Top points must be ordered before bottom points.')


def text_to_json(file, new_format=True):
    """
    Read in old or new ShakeMap 3 textfile rupture format and convert to
    GeoJSON.

    This will handle ShakeMap3.5-style fault text files, which can have the
    following format:
     - # at the top indicates a reference.
     - Lines beginning with a > indicate the end of one segment and the
       beginning of another.
     - Coordinates are specified in lat,lon,depth order.
     - Coordinates can be separated by commas or spaces.
     - Vertices can be specified in top-edge or bottom-edge first order.

    Args:
        file (str):
            Path to rupture file OR file-like object in GMT
            psxy format, where:

                * Rupture vertices are space/comma separated lat, lon, depth
                  triplets on a single line.
                * Rupture groups are separated by lines containing ">"
                * Rupture groups must be closed.
                * Verticies within a rupture group must start along the top
                  edge and move in the strike direction then move to the bottom
                  edge and move back in the opposite direction.

        new_format (bool):
            Indicates whether text rupture format is
            "old" (lat, lon, depth) or "new" (lon, lat, depth) style.

    Returns:
        dict: GeoJSON rupture dictionary.

    """
    isfile = False
    if hasattr(file, 'read'):
        f = file
    else:
        f = open(file, 'rt', encoding="latin-1")
        isfile = True

    reference = ''
    polygons = []
    polygon = []
    for line in f.readlines():
        if not len(line.strip()):
            continue

        if line.strip().startswith('#'):
            # Get reference string
            reference += line.strip().replace('#', '')
            continue

        if line.strip().startswith('>'):
            if not len(polygon):
                continue
            polygons.append(polygon)
            polygon = []
            continue

        # first try to split on whitespace
        parts = line.split()
        if len(parts) == 1:
            if new_format:
                raise ShakeLibException(
                    'Rupture file %s has unspecified delimiters.' % file)
            parts = line.split(',')
            if len(parts) == 1:
                raise ShakeLibException(
                    'Rupture file %s has unspecified delimiters.' % file)

        if len(parts) != 3:
            msg = 'Rupture file %s is not in lat, lon, depth format.'
            if new_format:
                'Rupture file %s is not in lon, lat, depth format.'
            raise ShakeLibException(msg % file)

        parts = [float(p) for p in parts]
        if not new_format:
            old_parts = parts.copy()
            parts[0] = old_parts[1]
            parts[1] = old_parts[0]
        polygon.append(parts)

    if len(polygon):
        polygons.append(polygon)

    if isfile:
        f.close()

    # Try to fix polygons
    original_polygons = polygons.copy()
    fixed = []
    n_polygons = len(polygons)
    for i in range(n_polygons):
        n_verts = len(polygons[i])
        success = False
        for j in range(n_verts - 1):
            try:
                _check_polygon(polygons[i])
                success = True
                break
            except ValueError:
                polygons[i] = _rotate_polygon(polygons[i])
        if success:
            fixed.append(True)
        else:
            fixed.append(False)

    if not all(fixed):
        polygons = original_polygons

    json_dict = {
        "type": "FeatureCollection",
        "metadata": {
            'reference': reference
        },
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "rupture type": "rupture extent"
                },
                "geometry": {
                    "type": "MultiPolygon",
                    "coordinates": [polygons]
                }
            }
        ]
    }
    validate_json(json_dict)

    return json_dict


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
            _check_polygon(polygons[i])


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
    f = d['features'][0]
    geom = f['geometry']
    if geom['type'] == 'Point':
        return False
    polygons = geom['coordinates'][0]
    try:
        len(polygons)
    except Exception:
        return False
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
                return False

            # Within each quad, top and bottom edges must be horizontal
            tops = np.array([quad[0].depth, quad[1].depth])
            if not np.isclose(tops[0], tops, rtol=0,
                              atol=constants.DEPTH_TOL).all():
                return False
            bots = np.array([quad[2].depth, quad[3].depth])
            if not np.isclose(bots[0], bots, rtol=0,
                              atol=constants.DEPTH_TOL).all():
                return False

    return True
