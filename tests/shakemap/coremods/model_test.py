#!/usr/bin/env python

import os
import os.path

import numpy as np

from shakemap.utils.config import get_config_paths
from shakemap.coremods.model import ModelModule
from shakemap.coremods.assemble import AssembleModule
from shakemap.coremods.plotregr import PlotRegr
from impactutils.io.smcontainers import ShakeMapOutputContainer
from common import clear_files, set_files

########################################################################
# Test sm_model
########################################################################


def test_model_2():

    #
    # This is a small grid with station data and dyfi data (should succeed)
    #
    install_path, data_path = get_config_paths()
    event_path = os.path.join(data_path, 'nc72282711', 'current')
    set_files(event_path, {'event.xml': 'event.xml',
                           'stationlist.xml.small': 'stationlist.xml',
                           'dyfi_dat.xml.small': 'dyfi_dat.xml',
                           'model.conf': 'model.conf',
                           'boat_fault.txt': 'boat_fault.txt'})
    assemble = AssembleModule('nc72282711', comment='Test comment.')
    assemble.execute()
    del assemble
    model = ModelModule('nc72282711')
    model.execute()
    del model
    #
    # Since we've done this, we might as well run plotregr, too
    #
    plotregr = PlotRegr('nc72282711')
    plotregr.execute()
    plotregr.writeContents()
    del plotregr
    clear_files(event_path)


def test_model_3():

    #
    # This is a small grid with DYFI data only (should succeed)
    #
    install_path, data_path = get_config_paths()
    event_path = os.path.join(data_path, 'nc72282711', 'current')
    set_files(event_path, {'event.xml': 'event.xml',
                           'dyfi_dat.xml.small': 'dyfi_dat.xml',
                           'model.conf': 'model.conf',
                           'boat_fault.txt': 'boat_fault.txt'})
    assemble = AssembleModule('nc72282711', comment='Test comment.')
    assemble.execute()
    model = ModelModule('nc72282711')
    model.execute()
    clear_files(event_path)


def test_model_4():

    #
    # Run with no data and no fault, and use the default extent.
    #
    install_path, data_path = get_config_paths()
    event_path = os.path.join(data_path, 'nc72282711', 'current')
    set_files(event_path, {
        'event.xml': 'event.xml',
        'model.conf': 'model.conf'
    })
    assemble = AssembleModule('nc72282711',
                              comment='Test comment.')
    assemble.execute()
    model = ModelModule('nc72282711')
    model.execute()
    clear_files(event_path)


def test_model_sim():

    #
    # Run with no data and no fault, and use the default extent.
    #
    install_path, data_path = get_config_paths()
    # event_path = os.path.join(data_path, 'planet9', 'current')
    assemble = AssembleModule('planet9',
                              comment='Test comment.')
    assemble.execute()
    model = ModelModule('planet9')
    model.execute()


def test_directivity():

    #
    # Turned on directivity in model.conf
    #
    install_path, data_path = get_config_paths()
    event_path = os.path.join(data_path, 'directivity_test', 'current')
    set_files(event_path, {
        'event.xml': 'event.xml',
        'model.conf': 'model.conf',
        'dir_fault.txt': 'dir_fault.txt'
    })
    assemble = AssembleModule('directivity_test',
                              comment='Test comment.')
    assemble.execute()
    model = ModelModule('directivity_test')
    model.execute()
    clear_files(event_path)
    hdf_file = os.path.join(event_path, 'products', 'shake_result.hdf')
    oc = ShakeMapOutputContainer.load(hdf_file)
    sa3 = np.exp(oc.getIMTGrids(
        'SA(3.0)', 'GREATER_OF_TWO_HORIZONTAL')['mean'])
    # np.testing.assert_allclose(np.max(sa3), 1.15864273)
    np.testing.assert_allclose(np.max(sa3), 1.1567265149442174)
    # np.testing.assert_allclose(np.min(sa3), 0.9278920)
    np.testing.assert_allclose(np.min(sa3), 0.88508818541678)
    oc.close()


def test_masking():
    install_path, data_path = get_config_paths()
    event_path = os.path.join(data_path, 'masking_test', 'current')
    set_files(event_path, {
        'event.xml': 'event.xml',
        'model.conf': 'model.conf',
        'au_continental_shelf.geojson': 'au_continental_shelf.geojson',
    })
    assemble = AssembleModule('masking_test',
                              comment='Test comment.')
    assemble.execute()
    model = ModelModule('masking_test')
    model.execute()
    clear_files(event_path)
    hdf_file = os.path.join(event_path, 'products', 'shake_result.hdf')
    oc = ShakeMapOutputContainer.load(hdf_file)
    sa3 = oc.getIMTGrids('SA(3.0)', 'GREATER_OF_TWO_HORIZONTAL')['mean']
    removed = np.isnan(sa3).astype(int)
    assert(removed[240, 240] == 1)
    assert(removed[260, 240] == 0)
    np.testing.assert_equal(removed[::100, ::100], [
        [1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 0, 1],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0]
    ])
    oc.close()


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_model_2()
    test_model_3()
    test_model_4()
    test_model_sim()
    test_directivity()
    test_masking()
