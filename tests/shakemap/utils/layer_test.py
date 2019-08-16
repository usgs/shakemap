#!/usr/bin/env python

import os.path
import pytest

from configobj import ConfigObj
import numpy as np
from shapely.geometry import Point

from shakemap.utils.layers import (get_layer_distances,
                                   dist_to_layer,
                                   update_config_regions,
                                   validate_config)

from shakemap.utils.config import (get_data_path,
                                   get_config_paths)

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))


def layers_equal(layer1, layer2):
    assert sorted(layer1.keys()) == sorted(layer2.keys())
    assert np.allclose(sorted(layer1.values()), sorted(layer2.values()))


def test_layers():
    data_path = get_data_path()
    layer_path = os.path.join(data_path, 'layers')

    elon = -117.0
    elat = 33.0
    layer_distances = get_layer_distances(elon, elat, layer_path)
    reference = {'italy': 8710.04538291321, 'hawaii': 4009.810418951339,
                 'chile': 4803.501119602294, 'japan': 7462.090628325871,
                 'europe_share': 6337.760076297577,
                 'induced': 1581.7186730857293,
                 'australia': 9822.519908573386, 'turkey': 9595.056080385448,
                 'greece': 9513.190936814879, 'china': 8749.15702828136,
                 'california': 0.0, 'new_zealand': 8859.317797859405,
                 'taiwan': 9698.206101592557}
    layers_equal(layer_distances, reference)

    elon = -97.5
    elat = 36.5
    layer_distances = get_layer_distances(elon, elat, layer_path)
    reference = {'italy': 7600.631037485594, 'hawaii': 5619.745326956643,
                 'chile': 3585.560465124193, 'japan': 8221.294305532281,
                 'europe_share': 5178.578550706942, 'induced': 0.0,
                 'australia': 10877.954598395638, 'turkey': 8706.788113517472,
                 'greece': 8551.021519240006, 'china': 9045.63060410134,
                 'california': 1511.4662456781018,
                 'new_zealand': 9950.334279593713,
                 'taiwan': 10301.898766277141}
    layers_equal(layer_distances, reference)

    elon = 121.0
    elat = 22.5
    layer_distances = get_layer_distances(elon, elat, layer_path)
    reference = {'italy': 8555.0797075834, 'hawaii': 7467.418933530909,
                 'chile': 12071.522116104166, 'japan': 1230.4098709858458,
                 'europe_share': 6828.880422152432,
                 'induced': 10327.489296208103, 'australia': 3857.7430367906,
                 'turkey': 6892.820456944929, 'greece': 8039.666277993378,
                 'china': 284.61969009227863, 'california': 9064.755577907308,
                 'new_zealand': 7707.863992372902, 'taiwan': 0.0}
    layers_equal(layer_distances, reference)

    #
    # Test for geometry type exception in dist_to_layer by
    # handing it a Point rather than a Polygon or MultiPolygon
    #
    p = Point()
    with pytest.raises(TypeError):
        dist_to_layer(0.0, 0.0, p)

    #
    # Test the updates to the config based on being in a layer (or not)
    #
    install_path, data_path = get_config_paths()
    global_data_path = os.path.join(os.path.expanduser('~'), 'shakemap_data')
    config = ConfigObj(os.path.join(install_path, 'config', 'select.conf'))
    validate_config(config, install_path, data_path, global_data_path)
    # Taiwan
    elon = 121.0
    elat = 22.5

    config = update_config_regions(elat, elon, config)
    assert config['tectonic_regions']['acr']['gmpe'] == \
        ['active_crustal_taiwan', 'active_crustal_taiwan_deep']

    config = ConfigObj(os.path.join(install_path, 'config', 'select.conf'))
    validate_config(config, install_path, data_path, global_data_path)
    # Induced
    elon = -97.5
    elat = 36.5

    config = update_config_regions(elat, elon, config)
    assert config['tectonic_regions']['scr']['gmpe'] == \
        ['stable_continental_nshmp2014_rlme', 'stable_continental_deep']

    config = ConfigObj(os.path.join(install_path, 'config', 'select.conf'))
    validate_config(config, install_path, data_path, global_data_path)
    # Not in a layer
    elon = -77.5
    elat = 36.5

    config = update_config_regions(elat, elon, config)
    assert config['tectonic_regions']['acr']['gmpe'] == \
        ['active_crustal_nshmp2014', 'active_crustal_deep']
    assert config['tectonic_regions']['scr']['gmpe'] == \
        ['stable_continental_nshmp2014_rlme', 'stable_continental_deep']

# def test_get_probability():
#     x1 = 0.0
#     p1 = 0.9
#     x2 = 20.0
#     p2 = 0.1

#     # test maximum probability
#     y1 = get_probability(0.0,x1,p1,x2,p2)
#     assert y1 == 0.9

#     # test minimum probability
#     y2 = get_probability(20.0,x1,p1,x2,p2)
#     assert y2 == 0.1

#     # make sure that minimum probability is a floor
#     y3 = get_probability(40.0,x1,p1,x2,p2)
#     assert y3 == 0.1

#     # test probability function
#     y4 = get_probability(10.0,x1,p1,x2,p2)
#     np.testing.assert_almost_equal(y4,0.5)

# def test_get_subduction_probabilities():
#     # we don't have slab depth or kagan angle, but we're crustal
#     results = {'SlabModelDepth':np.nan,
#                'TectonicSubtype':'ACR',
#                'KaganAngle':np.nan,
#                'SlabModelDepthUncertainty':np.nan,
#                'TensorType':'composite',
#                'FocalMechanism':'ALL'}
#     crustal,interface,intraslab = get_subduction_probabilities(results,0.0)
#     assert crustal == 0.8
#     assert interface == 0.1
#     assert intraslab == 0.1

#     # we don't have slab depth or kagan angle, but we're interface
#     results['TectonicSubtype'] = 'SZInter'
#     crustal,interface,intraslab = get_subduction_probabilities(results,20.0)
#     assert crustal == 0.1
#     assert interface == 0.8
#     assert intraslab == 0.1

#     # we don't have slab depth or kagan angle, but we're intraslab
#     results['TectonicSubtype'] = 'SZIntra'
#     crustal,interface,intraslab = get_subduction_probabilities(results,100.0)
#     assert crustal == 0.1
#     assert interface == 0.1
#     assert intraslab == 0.8

#     # we do have slab depth, no kagan angle, and RS mechanism
#     results = {'SlabModelDepth':20.0,
#                'TectonicSubtype':'ACR',
#                'KaganAngle':np.nan,
#                'SlabModelDepthUncertainty':10.0,
#                'TensorType':'composite',
#                'FocalMechanism':'RS'}
#     crustal,interface,intraslab = get_subduction_probabilities(results,20.0)
#     np.testing.assert_almost_equal(crustal,0.1625)
#     np.testing.assert_almost_equal(interface,0.675)
#     np.testing.assert_almost_equal(intraslab,0.1625)

#     # we do have slab depth, no kagan angle, and non-RS mechanism
#     results = {'SlabModelDepth':20.0,
#                'TectonicSubtype':'ACR',
#                'KaganAngle':np.nan,
#                'SlabModelDepthUncertainty':10.0,
#                'TensorType':'composite',
#                'FocalMechanism':'ALL'}
#     crustal,interface,intraslab = get_subduction_probabilities(results,20.0)
#     np.testing.assert_almost_equal(crustal,0.275)
#     np.testing.assert_almost_equal(interface,0.45)
#     np.testing.assert_almost_equal(intraslab,0.275)

#     # we have slab depth and kagan angle
#     results = {'SlabModelDepth':20.0,
#                'TectonicSubtype':'ACR',
#                'KaganAngle':0.0,
#                'SlabModelDepthUncertainty':10.0,
#                'TensorType':'composite',
#                'FocalMechanism':'ALL'}
#     crustal,interface,intraslab = get_subduction_probabilities(results,20.0)
#     np.testing.assert_almost_equal(crustal,0.05)
#     np.testing.assert_almost_equal(interface,0.9)
#     np.testing.assert_almost_equal(intraslab,0.05)

#     # we have slab depth and kagan angle, not right on the slab
#     results = {'SlabModelDepth':20.0,
#                'TectonicSubtype':'ACR',
#                'KaganAngle':0.0,
#                'SlabModelDepthUncertainty':10.0,
#                'TensorType':'composite',
#                'FocalMechanism':'ALL'}
#     crustal,interface,intraslab = get_subduction_probabilities(results,10.0)
#     np.testing.assert_almost_equal(crustal,0.25)
#     np.testing.assert_almost_equal(interface,0.5)
#     np.testing.assert_almost_equal(intraslab,0.25)

#     # we have slab depth and kagan angle, on the slab but kagan = 30
#     results = {'SlabModelDepth':20.0,
#                'TectonicSubtype':'ACR',
#                'KaganAngle':30.0,
#                'SlabModelDepthUncertainty':10.0,
#                'TensorType':'composite',
#                'FocalMechanism':'ALL'}
#     crustal,interface,intraslab = get_subduction_probabilities(results,20.0)
#     np.testing.assert_almost_equal(crustal,0.2525)
#     np.testing.assert_almost_equal(interface,0.495)
#     np.testing.assert_almost_equal(intraslab,0.2525)

#     # we have slab depth and kagan angle, off the slab and kagan = 30
#     results = {'SlabModelDepth':20.0,
#                'TectonicSubtype':'ACR',
#                'KaganAngle':30.0,
#                'SlabModelDepthUncertainty':10.0,
#                'TensorType':'composite',
#                'FocalMechanism':'ALL'}
#     crustal,interface,intraslab = get_subduction_probabilities(results,10.0)
#     np.testing.assert_almost_equal(crustal,0.3625,decimal=4)
#     np.testing.assert_almost_equal(interface,0.275,decimal=4)
#     np.testing.assert_almost_equal(intraslab,0.3625,decimal=4)


if __name__ == '__main__':
    test_layers()
