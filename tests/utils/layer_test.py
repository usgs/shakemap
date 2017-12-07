#!/usr/bin/env python

import os.path
import pytest

import numpy as np
from shapely.geometry import Point

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))

from shakemap.utils.layers import get_layer_distances, dist_to_layer
from shakemap.utils.config import get_data_path


def layers_equal(layer1, layer2):
    assert sorted(layer1.keys()) == sorted(layer2.keys())
    assert np.allclose(sorted(layer1.values()), sorted(layer2.values()))


def test_layers():
    data_path = get_data_path()
    layer_path = os.path.join(data_path, 'layers')

    elon = -117.0
    elat = 33.0
    layer_distances = get_layer_distances(elon, elat, layer_path)
    reference = {'induced': 1578.3879076203307,
                 'japan': 7972.1138613743387,
                 'taiwan': 11022.339157753582,
                 'california': 0.0}
    layers_equal(layer_distances, reference)

    elon = -97.5
    elat = 36.5
    layer_distances = get_layer_distances(elon, elat, layer_path)
    reference = {'induced': 0.0,
                 'japan': 8935.9779110700729,
                 'taiwan': 11997.837464370788,
                 'california': 1508.2155746648657}
    layers_equal(layer_distances, reference)

    elon = 121.0
    elat = 22.5
    layer_distances = get_layer_distances(elon, elat, layer_path)
    reference = {'induced': 12041.424518656486,
                 'japan': 1231.8954391427453,
                 'taiwan': 0.0,
                 'california': 10085.281293655946}
    layers_equal(layer_distances, reference)

    #
    # Test for geometry type exception in dist_to_layer by
    # handing it a Point rather than a Polygon or MultiPolygon
    #
    p = Point()
    with pytest.raises(TypeError):
        distance = dist_to_layer(0.0, 0.0, p)


if __name__ == '__main__':
    test_layers()
