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

    datadir = os.path.join(
        homedir, '..', '..', 'data', 'eventdata', 'us1000db40', 'current')
    eventxml = os.path.join(datadir, 'event.xml')
#    momentfile = os.path.join(datadir, 'moment.xml')
    org = Origin.fromFile(eventxml)
    # It seems the following line does not work and needs to get fixed in
    # select.py so that the moment tensor gets used.
    # tensor_params = org.moment.copy()
    gmpe_list, weight_list, strec_results = get_weights(
        org, config, tensor_params=None)
    assert len(gmpe_list) == 3
    assert len(weight_list) == 3
    assert gmpe_list[0] == 'subduction_crustal'
    assert gmpe_list[1] == 'subduction_interface_nshmp2014'
    assert gmpe_list[2] == 'subduction_slab_nshmp2014'
    np.testing.assert_allclose(weight_list[0], 0.00543576, rtol=1e-05)
    np.testing.assert_allclose(weight_list[1], 0.96074659, rtol=1e-05)
    np.testing.assert_allclose(weight_list[2], 0.03381765, rtol=1e-05)


if __name__ == '__main__':
    test_probs()
