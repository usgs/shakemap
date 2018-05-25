#!/usr/bin/env python

import os
import os.path

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.model import ModelModule
from shakemap.coremods.assemble import AssembleModule
from shakemap.coremods.plotregr import PlotRegr

########################################################################
# Test sm_model
########################################################################


def test_model_1():

    installpath, datapath = get_config_paths()
    #
    # This is Northridge for a set of output points (not a grid)
    #
    assemble = AssembleModule('northridge_points', comment='Test comment.')
    assemble.execute()
    model = ModelModule('northridge_points')
    model.execute()


def test_model_2():

    #
    # This is a small grid with station data and dyfi data (should succeed)
    #
    assemble = AssembleModule('nc72282711', comment='Test comment.')
    assemble.execute()
    model = ModelModule('nc72282711')
    model.execute()
    #
    # Since we've done this, we might as well run plotregr, too
    #
    plotregr = PlotRegr('nc72282711')
    plotregr.execute()
    plotregr.writeContents()
    pass


def test_model_3():

    #
    # This is a small grid with DYFI data only (should succeed)
    #
    assemble = AssembleModule('nc72282711_dyfi', comment='Test comment.')
    assemble.execute()
    model = ModelModule('nc72282711_dyfi')
    model.execute()


def test_model_4():

    #
    # Run with no data and no fault, and use the default extent.
    #
    assemble = AssembleModule('nc72282711_nodata_nofault',
                              comment='Test comment.')
    assemble.execute()
    model = ModelModule('nc72282711_nodata_nofault')
    model.execute()


def test_model_5():

    #
    # Set the bias and outlier magnitude limits low to test additional
    # code branches
    #
    assemble = AssembleModule('nc72282711_nofault', comment='Test comment.')
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
