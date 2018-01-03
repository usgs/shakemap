import glob
import os.path
from functools import partial

# third party imports
import pyproj
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import transform
import shapely.wkt
import numpy as np
from openquake.hazardlib.geo.geodetic import min_distance_to_segment
from strec.subtype import SubductionSelector

def nearest_edge(elon, elat, poly):
    """
    Return the distance from a point to the nearest edge of a
    polygon.

    Args:
        elon (float): The longitude of the reference point.
        elat (float): The latitude of the reference point.
        poly (Polygon): An instance of a shapely Polygon.

    Returns:
        float: The distance (in km) from the reference point to the
        nearest edge (or vertex) of the polygon.
    """
    elon_arr = np.array([elon])
    elat_arr = np.array([elat])
    x, y = poly.exterior.xy
    nearest = 99999.
    for ix in range(1, len(x) - 1):
        dd = min_distance_to_segment(np.array(x[ix - 1:ix + 1]),
                                     np.array(y[ix - 1:ix + 1]),
                                     elon_arr, elat_arr)
        if np.abs(dd[0]) < nearest:
            nearest = dd[0]
    return nearest


def dist_to_layer(elon, elat, geom):
    """
    Return the distance from a point to the polygon(s) in a layer; zero if
    the point is inside the polygon. If the nearest edge of the polygon is
    greater than 5000 km from the point, the point cannot be inside the
    polygon and the distance reported will be the distance to the nearest
    edge. So don't make polygons too big.

    Args:
        elon (float): The longitude of the reference point.
        elat (float): The latitude of the reference point.
        geom (Polygon or MultiPolygon): An instance of a shapely Polygon
            or MultiPolygon.

    Returns:
        float: The distance (in km) from the reference point to the
        nearest polygon in the layer. The distance will be zero if the
        point lies inside the polygon.
    """
    if isinstance(geom, Polygon):
        plist = [geom]
    elif isinstance(geom, MultiPolygon):
        plist = list(geom.geoms)
    else:
        raise TypeError('Invalid geometry type in layer: %s' % type(geom))

    project = partial(
        pyproj.transform,
        pyproj.Proj(proj='latlong'),
        pyproj.Proj(proj='aeqd  +lat_0=%f +lon_0=%f +R=6371' % (elat, elon)))
    ep = Point(0.0, 0.0)
    min_dist = 99999.
    for poly in plist:
        nearest = nearest_edge(elon, elat, poly)
        if nearest < 5000:
            nearest = ep.distance(transform(project, poly))
        if nearest < min_dist:
            min_dist = nearest
        if min_dist == 0:
            break
    return min_dist


def get_layer_distances(elon, elat, layer_dir):
    """
    Return the distances from a point to the nearest polygon in each
    layer file found in 'layer_dir'. The distance will be zero if
    the point is inside a polygon. If the nearest edge of a polygon is
    greater than 5000 km from the point, the point cannot be inside the
    polygon and the distance reported will be the distance to the
    nearest edge. So don't make polygons too big.

    The layer files should be written in Well-Known Text (with .wkt
    extensions), and should contain either a single POLYGON or
    MULTIPOLYGON object. The layer name will be the file's basename.

    Args:
        elon (float): The longitude of the reference point.
        elat (float): The latitude of the reference point.
        layer_dir (str): The path to the directory containg the layer
            files.

    Returns:
        dict: A dictionary where the keys are the layer names, and the
        values are the distance (in km) from the reference point to the
        nearest polygon in the layer. The distance will be zero if the
        point lies inside the polygon.
    """
    layer_files = glob.glob(os.path.join(layer_dir, '*.wkt'))
    dist_dict = {}
    for file in layer_files:
        layer_name = os.path.splitext(os.path.basename(file))[0]
        with open(file, 'r') as fd:
            data = fd.read()
        geom = shapely.wkt.loads(data)
        dist_dict[layer_name] = dist_to_layer(elon, elat, geom)
    return dist_dict

def get_subduction_probabilities(results,depth):
    # inputs to algorithm
    kagan = results['KaganAngle'] #can be nan
    slab_depth = results['SlabModelDepth']
    slab_depth_error = results['SlabModelDepthUncertainty']
    focal = results['FocalMechanism']
    tensor_type = results['TensorType']
    subtype = results['TectonicSubtype']
    if np.isnan(kagan):
        if np.isnan(slab_depth):
            if subtype == 'SZInter':
                crustal_prob = 0.1
                intraslab_prob = 0.1
                interface_prob = 0.8
            elif subtype == 'ACR':
                crustal_prob = 0.8
                intraslab_prob = 0.1
                interface_prob = 0.1
            else:
                crustal_prob = 0.1
                intraslab_prob = 0.8
                interface_prob = 0.1
        else:
            # here we don't have kagan angle, but we do have depth to slab
            

def get_tectonic_regions(elon, elat, edepth, eid):
    selector = SubductionSelector()
    results = selector.getSubductionTypeByID(eid)
    strec_out = {}
    strec_out['focal_mech'] = results['FocalMechanism']

    #figure out the probabilities of subduction zone
    crustal_prob = 0.0
    interface_prob = 0.0
    intraslab_prob = 0.0
    if results['TectonicRegion'] == 'Subduction':
        crustal_prob,interface_prob,intraslab_prob = get_subduction_probabilities(results,depth)
    
    regions = {'acr':{'distance':results['DistanceToActive']},
               'scr':{'distance':results['DistanceToStable']},
               'volcanic':{'distance':results['DistanceToVolcanic']},
               'subduction':{'distance':results['DistanceToSubduction'],
                             'probabilities':{'crustal':crustal_prob,
                                              'interface':interface_prob,
                                              'intraslab':intraslab_prob}}}
    
    strec_out['tectonic_regions'] = regions
    # strec_out = {
    #     'focal_mech': 'ALL',
    #     'tectonic_regions': {
    #         'acr': {
    #             'distance': 0.0,
    #         },
    #         'scr': {
    #             'distance': 1500.0,
    #         },
    #         'subduction': {
    #             'distance': 1400.0,
    #             'probabilities': {
    #                 'crustal': 0,
    #                 'interface': 0,
    #                 'intraslab': 0,
    #             }
    #         },
    #         'volcanic': {
    #             'distance': 2300.0,
    #         },
    #     }
    # }

    return strec_out
