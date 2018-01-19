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

def get_probability(x,x1,p1,x2,p2):
    if x <= x1:
        prob = p1
    elif x >= x2:
        prob = p2
    else:
        slope = (p1-p2)/(x1-x2)
        intercept = p1 - slope * x1
        prob = x * slope + intercept
    return prob

def get_subduction_probabilities(results,depth,config):
    
    # inputs to algorithm from STREC

    # the angle between moment tensor and slab
    kagan = results['KaganAngle'] #can be nan

    # Depth to slab
    slab_depth = results['SlabModelDepth']

    # Error in depth to slab
    slab_depth_error = results['SlabModelDepthUncertainty']

    # what is the effective bottom of the interface zone?
    max_interface_depth = results['SlabModelMaximumDepth']

    # Calculate the probability of interface given the
    # (absolute value of) difference between hypocenter and depth to slab.
    dz = np.abs(depth - slab_depth)
    x1 = config['p_int_hypo']['x1'] + slab_depth_error
    x2 = config['p_int_hypo']['x2'] + slab_depth_error
    p1 = config['p_int_hypo']['p1']
    p2 = config['p_int_hypo']['p2']
    p_int_hypo = get_probability(dz,x1,p1,x2,p2)

    # Calculate probability of interface given Kagan's angle
    if np.isfinite(kagan):
        x1 = config['p_int_kagan']['x1']
        x2 = config['p_int_kagan']['x2']
        p1 = config['p_int_kagan']['p1']
        p2 = config['p_int_kagan']['p2']
        p_int_kagan = get_probability(kagan,x1,p1,x2,p2)
    else:
        p_int_kagan = config['p_kagan_default']

    # Calculate probability that event occurred above bottom of seismogenic
    # zone, given to us by the Slab model.
    x1 = max_interface_depth + config['p_int_sz']['x1']
    x2 = max_interface_depth + config['p_int_sz']['x2']
    p1 = config['p_int_sz']['p1']
    p2 = config['p_int_sz']['p2']
    p_int_sz = get_probability(depth,x1,p1,x2,p2)

    # Calculate combined probability of interface
    p_int = p_int_hypo * p_int_kagan * p_int_sz

    # Calculate probability that the earthquake lies above the slab
    # and is thus crustal.
    x1 = config['p_crust_slab']['x1']
    x2 = config['p_crust_slab']['x2']
    p1 = config['p_crust_slab']['p1']
    p2 = config['p_crust_slab']['p2']
    p_crust_slab = get_probability((depth-slab_depth),x1,p1,x2,p2)

    # Calculate probability that the earthquake lies within the crust
    x1 = config['p_crust_hypo']['x1']
    x2 = config['p_crust_hypo']['x2']
    p1 = config['p_crust_hypo']['p1']
    p2 = config['p_crust_hypo']['p2']
    p_crust_hypo = get_probability(depth,x1,p1,x2,p2)
    
    # Calculate probability of crustal
    p_crustal = (1 - p_int) * p_crust_slab * p_crust_hypo

    # Calculate probability of intraslab
    p_slab = 1 - (p_int + p_crustal)
        

    probs = {'crustal_probability' : p_crustal,
             'interface_probability' : p_int,
             'intraslab_probability' : p_slab,
             'depth_interface_probability' : p_int_hypo,
             'kagan_interface_probability' : p_int_kagan,
             'seismo_zone_interface_probability' : p_int_sz,
             'above_slab_probability' : p_crust_slab,
             'within_slab_probability' : p_crust_hypo}
    
    return probs

def get_tectonic_regions(elon, elat, edepth, eid,config):
    selector = SubductionSelector()
    results = selector.getSubductionTypeByID(eid)
    strec_out = {}
    strec_out['focal_mech'] = results['FocalMechanism']

    #figure out the probabilities of subduction zone
    crustal_prob = 0.0
    interface_prob = 0.0
    intraslab_prob = 0.0
    depth_interface_prob = 0.0
    kagan_interface_prob = 0.0
    sz_interface_prob = 0.0
    crust_slab_prob = 0.0
    within_slab_prob = 0.0
    if not np.isnan(results['SlabModelDepth']):
        subcfg = config['subduction']
        probs = get_subduction_probabilities(results,edepth,subcfg)
        crustal_prob = probs['crustal_probability']
        interface_prob = probs['interface_probability']
        intraslab_prob = probs['intraslab_probability']
        depth_interface_prob = probs['depth_interface_probability']
        kagan_interface_prob = probs['kagan_interface_probability']
        sz_interface_prob = probs['seismo_zone_interface_probability']
        crust_slab_prob = probs['above_slab_probability']
        within_slab_prob = probs['within_slab_probability']
        
         
    sub_probs = {'crustal':crustal_prob,
                 'interface':interface_prob,
                 'intraslab':intraslab_prob,
                 'depth_interface':depth_interface_prob,
                 }

    sub_inputs = {'depth' : edepth,
                  'slab_depth' : results['SlabModelDepth'],
                  'slab_depth_error' : results['SlabModelDepthUncertainty'] ,
                  'kagan_angle' : results['KaganAngle'],
                  'focal_mech' : results['FocalMechanism']}
        
    regions = {'acr':{'distance':results['DistanceToActive']},
               'scr':{'distance':results['DistanceToStable']},
               'volcanic':{'distance':results['DistanceToVolcanic']},
               'subduction':{'distance':results['DistanceToSubduction'],
                             'probabilities':sub_probs,
                             'inputs':sub_inputs}}
    
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
