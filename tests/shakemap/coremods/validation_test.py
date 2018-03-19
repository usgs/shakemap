#!/usr/bin/env python
import os
import os.path

import matplotlib
matplotlib.use('Agg')

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.model import ModelModule
from shakemap.coremods.assemble import AssembleModule
from shakemap.coremods.xtestplot import XTestPlot
from shakemap.coremods.xtestplot_spectra import XTestPlotSpectra
from shakemap.coremods.xtestimage import XTestImage
from shakemap.coremods.plotregr import PlotRegr

########################################################################
# Test the nullgmpe and the dummy correlation function as well as
# the two xtestplot functions
########################################################################


def test_validation():

    installpath, datapath = get_config_paths()
    try:
        #
        # Test xtestplot on validation event 0006
        #
        assemble = AssembleModule('validation_test_0006',
                                  comment='Test comment.')
        assemble.execute()
        model = ModelModule('validation_test_0006')
        model.execute()
        plot = XTestPlot('validation_test_0006')
        plot.execute()

        #
        # Test xtestplot_spectra on validation event 0007
        #
        assemble = AssembleModule('validation_test_0007',
                                  comment='Test comment.')
        assemble.execute()
        model = ModelModule('validation_test_0007')
        model.execute()
        plot = XTestPlotSpectra('validation_test_0007')
        plot.execute()

        #
        # Test xtestimage on validation event 0013
        #
        assemble = AssembleModule('validation_test_0013',
                                  comment='Test comment.')
        assemble.execute()
        model = ModelModule('validation_test_0013')
        model.execute()
        plot = XTestImage('validation_test_0013')
        plot.execute()
        regr = PlotRegr('validation_test_0013')
        regr.execute()
    finally:
        data_file = os.path.join(datapath, 'validation_test_0006', 'current',
                                 'shake_data.hdf')
        if os.path.isfile(data_file):
            os.remove(data_file)
        res_file = os.path.join(datapath, 'validation_test_0006', 'current',
                                'products', 'shake_results.hdf')
        if os.path.isfile(res_file):
            os.remove(res_file)
        data_file = os.path.join(datapath, 'validation_test_0007', 'current',
                                 'shake_data.hdf')
        if os.path.isfile(data_file):
            os.remove(data_file)
        res_file = os.path.join(datapath, 'validation_test_0007', 'current',
                                'products', 'shake_results.hdf')
        if os.path.isfile(res_file):
            os.remove(res_file)
        data_file = os.path.join(datapath, 'validation_test_0013', 'current',
                                 'shake_data.hdf')
        if os.path.isfile(data_file):
            os.remove(data_file)
        res_file = os.path.join(datapath, 'validation_test_0013', 'current',
                                'products', 'shake_results.hdf')
        if os.path.isfile(res_file):
            os.remove(res_file)


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_validation()
