#!/usr/bin/env python

import numpy as np
import os.path
from configobj import ConfigObj

from shakemap.utils.config import get_config_paths
from shakemap.utils.probs import get_weights
from shakemap.coremods.select import validate_config

from shakelib.rupture.origin import Origin

homedir = os.path.dirname(os.path.abspath(__file__))


def test_probs():
    install_path, config_path = get_config_paths()
    config = ConfigObj(os.path.join(install_path, 'config', 'select.conf'))
    validate_config(config, install_path)

    #
    # Real event: M 6.9 - 151km SSW of Kokopo, Papua New Guinea
    # 2018-03-29 21:25:36 UTC
    # Note: slab model depth is 33.4638 km
    #
    datadir = os.path.join(
        homedir, '..', '..', 'data', 'eventdata', 'us1000db40', 'current')
    eventxml = os.path.join(datadir, 'event.xml')
    org = Origin.fromFile(eventxml)
    gmpe_list, weight_list, strec_results = get_weights(org, config)
    assert len(gmpe_list) == 1
    assert len(weight_list) == 1
    assert gmpe_list[0] == 'subduction_interface_nshmp2014'
    np.testing.assert_allclose(weight_list[0], 1.0, rtol=1e-05)

    # move hypo to the surface:
    org.depth = 0
    gmpe_list, weight_list, strec_results = get_weights(org, config)
    # dz is 33.4638
    # slab uncertainty is 10, so x2 in int|abs(dz) is 29
    # so m_{int|abs(dz)} = 0.15
    np.testing.assert_allclose(
        weight_list,
        np.array([0.85,  0.15]),
        rtol=1e-05)
    assert gmpe_list[0] == 'subduction_crustal'
    assert gmpe_list[1] == 'subduction_interface_nshmp2014'

    # move hypo to be deeper:
    org.depth = 65
    gmpe_list, weight_list, strec_results = get_weights(org, config)
    # dz is 33.4638
    # slab uncertainty is 10, so x2 in int|abs(dz) is 29
    # so m_{int|abs(dz)} = 0.15
    np.testing.assert_allclose(
        weight_list,
        np.array([0.15,  0.85]),
        rtol=1e-05)
    assert gmpe_list[0] == 'subduction_interface_nshmp2014'
    assert gmpe_list[1] == 'subduction_slab_nshmp2014'


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_probs()
