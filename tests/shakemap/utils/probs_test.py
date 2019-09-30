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
    install_path, data_path = get_config_paths()
    global_data_path = os.path.join(os.path.expanduser('~'), 'shakemap_data')
    config = ConfigObj(os.path.join(install_path, 'config', 'select.conf'))
    validate_config(config, install_path, data_path, global_data_path)

    #
    # Real event: M 6.9 - 151km SSW of Kokopo, Papua New Guinea
    # 2018-03-29 21:25:36 UTC
    # Note: slab model depth is 33.4638 km
    #
    datadir = os.path.join(
        homedir, '..', '..', 'data', 'eventdata', 'us1000db40', 'current')
    eventxml = os.path.join(datadir, 'event.xml')
    org = Origin.fromFile(eventxml)
    gmmodel, strec_results = get_weights(org, config)
    assert len(gmmodel['gmpelist']) == 2
    assert len(gmmodel['weightlist']) == 2
    assert gmmodel['gmpelist'][0] == 'active_crustal_deep'
    assert gmmodel['gmpelist'][1] == 'subduction_interface_nshmp2014'
    np.testing.assert_allclose(gmmodel['weightlist'], [0.191049, 0.808951],
                               rtol=1e-05)
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    # move hypo to the surface:
    org.depth = 0
    gmmodel, strec_results = get_weights(org, config)
    # dz is 33.4638
    # slab uncertainty is 10, so x2 in int|abs(dz) is 29
    # so m_{int|abs(dz)} = 0.15
    np.testing.assert_allclose(
        gmmodel['weightlist'],
        np.array([0.191049, 0.32369088, 0.48526013]),
        rtol=1e-05)
    assert gmmodel['gmpelist'][0] == 'active_crustal_nshmp2014'
    assert gmmodel['gmpelist'][1] == 'subduction_crustal'
    assert gmmodel['gmpelist'][2] == 'subduction_interface_nshmp2014'
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    # move hypo to be deeper:
    org.depth = 65
    gmmodel, strec_results = get_weights(org, config)
    # dz is 33.4638
    # slab uncertainty is 10, so x2 in int|abs(dz) is 29
    # so m_{int|abs(dz)} = 0.15
    np.testing.assert_allclose(
        gmmodel['weightlist'],
        np.array([0.191049, 0.808951]), rtol=1e-05)
    assert gmmodel['gmpelist'][0] == 'active_crustal_deep'
    assert gmmodel['gmpelist'][1] == 'subduction_slab_nshmp2014'
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    #
    # Fake event -- subduction, but not within the Slab model
    #
    org.lon = -80.0
    org.lat = 9.0
    org.depth = 0
    gmmodel, strec_results = get_weights(org, config)
    assert gmmodel['weightlist'] == [1.0]
    assert gmmodel['gmpelist'] == ['subduction_crustal']
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    org.depth = 70.0
    gmmodel, strec_results = get_weights(org, config)
    assert gmmodel['weightlist'] == [1.0]
    assert gmmodel['gmpelist'] == ['subduction_slab_nshmp2014']
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    org.depth = 30.0
    gmmodel, strec_results = get_weights(org, config)
    assert gmmodel['weightlist'] == [1.0]
    assert gmmodel['gmpelist'] == ['subduction_interface_nshmp2014']
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    org.depth = 50.0
    gmmodel, strec_results = get_weights(org, config)
    assert gmmodel['weightlist'][0] == 0.5 and gmmodel['weightlist'][1] == 0.5
    assert gmmodel['gmpelist'] == ['subduction_interface_nshmp2014',
                                   'subduction_slab_nshmp2014']
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    org.depth = 20.0
    gmmodel, strec_results = get_weights(org, config)
    assert np.allclose(gmmodel['weightlist'], [0.7, 0.3])
    assert gmmodel['gmpelist'] == ['subduction_crustal',
                                   'subduction_interface_nshmp2014']
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    # Adding some new events that Bruce identified as being problematic
    datadir = os.path.join(
        homedir, '..', '..', 'data', 'eventdata', 'us1000ehd9', 'current')
    eventxml = os.path.join(datadir, 'event.xml')
    org = Origin.fromFile(eventxml)
    gmmodel, strec_results = get_weights(org, config)
    assert len(gmmodel['gmpelist']) == 4
    assert len(gmmodel['weightlist']) == 4
    np.testing.assert_allclose(
        np.sum(gmmodel['weightlist']), 1.0, rtol=1e-05)
    assert gmmodel['gmpelist'][0] == 'active_crustal_nshmp2014'
    assert gmmodel['gmpelist'][1] == 'active_crustal_deep'
    assert gmmodel['gmpelist'][2] == 'subduction_interface_nshmp2014'
    assert gmmodel['gmpelist'][3] == 'subduction_slab_nshmp2014'
    weights = np.array([0.14002268, 0.31823337, 0.13543599, 0.40630796])
    np.testing.assert_allclose(weights, gmmodel['weightlist'], rtol=1e-05)
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    datadir = os.path.join(
        homedir, '..', '..', 'data', 'eventdata', 'nc73027396', 'current')
    eventxml = os.path.join(datadir, 'event.xml')
    org = Origin.fromFile(eventxml)
    gmmodel, strec_results = get_weights(org, config)
    assert len(gmmodel['gmpelist']) == 2
    assert len(gmmodel['weightlist']) == 2
    np.testing.assert_allclose(
        np.sum(gmmodel['weightlist']), 1.0, rtol=1e-05)
    assert gmmodel['gmpelist'][0] == 'active_crustal_nshmp2014'
    assert gmmodel['gmpelist'][1] == 'subduction_crustal'
    weights = np.array([0.55958735, 0.44041265])
    np.testing.assert_allclose(weights, gmmodel['weightlist'], rtol=1e-05)
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    datadir = os.path.join(
        homedir, '..', '..', 'data', 'eventdata', 'nc73027746', 'current')
    eventxml = os.path.join(datadir, 'event.xml')
    org = Origin.fromFile(eventxml)
    gmmodel, strec_results = get_weights(org, config)
    assert len(gmmodel['gmpelist']) == 2
    assert len(gmmodel['weightlist']) == 2
    np.testing.assert_allclose(
        np.sum(gmmodel['weightlist']), 1.0, rtol=1e-05)
    assert gmmodel['gmpelist'][0] == 'active_crustal_nshmp2014'
    assert gmmodel['gmpelist'][1] == 'subduction_crustal'
    weights = np.array([0.65820838, 0.34179162])
    np.testing.assert_allclose(weights, gmmodel['weightlist'], rtol=1e-05)
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    datadir = os.path.join(
        homedir, '..', '..', 'data', 'eventdata', 'us1000er8i', 'current')
    eventxml = os.path.join(datadir, 'event.xml')
    org = Origin.fromFile(eventxml)
    gmmodel, strec_results = get_weights(org, config)
    assert len(gmmodel['gmpelist']) == 2
    assert len(gmmodel['weightlist']) == 2
    np.testing.assert_allclose(
        np.sum(gmmodel['weightlist']), 1.0, rtol=1e-05)
    assert gmmodel['gmpelist'][0] == 'active_crustal_nshmp2014'
    assert gmmodel['gmpelist'][1] == 'subduction_crustal'
    weights = np.array([0.59531667, 0.40468333])
    np.testing.assert_allclose(weights, gmmodel['weightlist'], rtol=1e-05)
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    datadir = os.path.join(
        homedir, '..', '..', 'data', 'eventdata', 'us1000e0j1', 'current')
    eventxml = os.path.join(datadir, 'event.xml')
    org = Origin.fromFile(eventxml)
    gmmodel, strec_results = get_weights(org, config)
    assert len(gmmodel['gmpelist']) == 2
    assert len(gmmodel['weightlist']) == 2
    np.testing.assert_allclose(
        np.sum(gmmodel['weightlist']), 1.0, rtol=1e-05)
    assert gmmodel['gmpelist'][0] == 'active_crustal_nshmp2014'
    assert gmmodel['gmpelist'][1] == 'subduction_crustal'
    weights = np.array([0.87367475, 0.12632525])
    np.testing.assert_allclose(weights, gmmodel['weightlist'], rtol=1e-05)
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'

    datadir = os.path.join(
        homedir, '..', '..', 'data', 'eventdata', 'us1000db5t', 'current')
    eventxml = os.path.join(datadir, 'event.xml')
    org = Origin.fromFile(eventxml)
    gmmodel, strec_results = get_weights(org, config)
    assert len(gmmodel['gmpelist']) == 2
    assert len(gmmodel['weightlist']) == 2
    np.testing.assert_allclose(
        np.sum(gmmodel['weightlist']), 1.0, rtol=1e-05)
    assert gmmodel['gmpelist'][0] == 'active_crustal_nshmp2014'
    assert gmmodel['gmpelist'][1] == 'subduction_crustal'
    weights = np.array([0.88676164, 0.11323836])
    np.testing.assert_allclose(weights, gmmodel['weightlist'], rtol=1e-05)
    assert gmmodel['ipe'] == 'VirtualIPE'
    assert gmmodel['gmice'] == 'WGRW12'
    assert gmmodel['ccf'] == 'LB13'


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_probs()
