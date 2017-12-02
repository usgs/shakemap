#!/usr/bin/env python

import os
import os.path
import shutil

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.model import ModelModule
from shakemap.coremods.assemble import AssembleModule

########################################################################
# Test sm_model
########################################################################
def test_model_1():

    installpath, datapath = get_config_paths()
    #
    # This is Northridge for a set of output points (not a grid)
    # Remove the products directory to hit the code that makes it
    # (should succeed)
    #
    assemble = AssembleModule('northridge_points')
    assemble.execute()
    products_dir = os.path.join(datapath, 'northridge_points', 'current', 
                            'products')
    if os.path.isdir(products_dir):
        shutil.rmtree(products_dir)
    model = ModelModule('northridge_points')
    model.execute()

def test_model_2():

    #
    # This is a small grid with station data only (should succeed)
    #
    assemble = AssembleModule('nc72282711')
    assemble.execute()
    model = ModelModule('nc72282711')
    model.execute()

def test_model_3():

    #
    # This is a small grid with DYFI data only (should succeed)
    #
    assemble = AssembleModule('nc72282711_dyfi')
    assemble.execute()
    model = ModelModule('nc72282711_dyfi')
    model.execute()

def test_model_4():

    #
    # Run with no data and no fault, and use the default extent.
    #
    assemble = AssembleModule('nc72282711_nodata_nofault')
    assemble.execute()
    model = ModelModule('nc72282711_nodata_nofault')
    model.execute()

def test_model_5():

    #
    # Set the bias and outlier magnitude limits low to test additional
    # code branches
    #
    assemble = AssembleModule('nc72282711_nofault')
    assemble.execute()
    model = ModelModule('nc72282711_nofault')
    model.execute()

def test_model_6():

    installpath, datapath = get_config_paths()
    #
    # This event exists, but we hide the input hdf file (should fail)
    #
    hdf_file = os.path.join(datapath, 'nc72282711_dyfi', 'current', 
                            'shake_data.hdf')
    if os.path.isfile(hdf_file):
        os.remove(hdf_file)
    model = ModelModule('nc72282711_dyfi')
    with pytest.raises(FileNotFoundError):
        model.execute()

def test_model_7():

    #
    # This event doesn't exist (should fail)
    #
    model = ModelModule('not_an_event')
    with pytest.raises(NotADirectoryError):
        model.execute()

if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_model_1()
    test_model_2()
    test_model_3()
    test_model_4()
    test_model_5()
    test_model_6()
    test_model_7()
