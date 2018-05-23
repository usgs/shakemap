#!/usr/bin/env python
import os
import os.path

from shakemap.utils.config import get_config_paths
from shakemap.coremods.model import ModelModule
from shakemap.coremods.assemble import AssembleModule
from shakemap.coremods.xtestplot import XTestPlot
from shakemap.coremods.xtestplot_spectra import XTestPlotSpectra
from shakemap.coremods.xtestimage import XTestImage
from shakemap.coremods.xtestplot_multi import XTestPlotMulti
from shakemap.coremods.plotregr import PlotRegr

########################################################################
# Test the nullgmpe and the dummy correlation function as well as
# the two xtestplot functions
########################################################################


def test_verification():

    installpath, datapath = get_config_paths()
    try:
        #
        # Test xtestplot on verification event 0006
        #
        assemble = AssembleModule('verification_test_0006',
                                  comment='Test comment.')
        assemble.execute()
        model = ModelModule('verification_test_0006')
        model.execute()
        plot = XTestPlot('verification_test_0006')
        plot.execute()

        #
        # Test xtestplot_spectra on verification event 0007
        #
        assemble = AssembleModule('verification_test_0007',
                                  comment='Test comment.')
        assemble.execute()
        model = ModelModule('verification_test_0007')
        model.execute()
        plot = XTestPlotSpectra('verification_test_0007')
        plot.execute()

        #
        # Test xtestimage on verification event 0011
        #
        assemble = AssembleModule('verification_test_0011',
                                  comment='Test comment.')
        assemble.execute()
        model = ModelModule('verification_test_0011')
        model.execute()
        plot = XTestImage('verification_test_0011')
        plot.execute()
        regr = PlotRegr('verification_test_0011')
        regr.execute()

        #
        # Test xtestplot_multi on event 0008x
        #
        for vt in ('8a', '8b', '8c', '8d', '8e'):
            assemble = AssembleModule('verification_test_000%s' % vt,
                                      comment='Test comment.')
            assemble.execute()
            model = ModelModule('verification_test_000%s' % vt)
            model.execute()
        plot = XTestPlotMulti('verification_test_0008')
        plot.execute()

    finally:
        data_file = os.path.join(datapath, 'verification_test_0006', 'current',
                                 'shake_data.hdf')
        if os.path.isfile(data_file):
            os.remove(data_file)
        res_file = os.path.join(datapath, 'verification_test_0006', 'current',
                                'products', 'shake_results.hdf')
        if os.path.isfile(res_file):
            os.remove(res_file)
        data_file = os.path.join(datapath, 'verification_test_0007', 'current',
                                 'shake_data.hdf')
        if os.path.isfile(data_file):
            os.remove(data_file)
        res_file = os.path.join(datapath, 'verification_test_0007', 'current',
                                'products', 'shake_results.hdf')
        if os.path.isfile(res_file):
            os.remove(res_file)
        data_file = os.path.join(datapath, 'verification_test_0011', 'current',
                                 'shake_data.hdf')
        if os.path.isfile(data_file):
            os.remove(data_file)
        res_file = os.path.join(datapath, 'verification_test_0011', 'current',
                                'products', 'shake_results.hdf')
        if os.path.isfile(res_file):
            os.remove(res_file)
        for vt in ('8a', '8b', '8c', '8d', '8e'):
            evid = 'verification_test_000%s' % vt
            data_file = os.path.join(datapath, evid, 'current',
                                     'shake_data.hdf')
            if os.path.isfile(data_file):
                os.remove(data_file)
            res_file = os.path.join(datapath, evid, 'current',
                                    'products', 'shake_results.hdf')
            if os.path.isfile(res_file):
                os.remove(res_file)


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_verification()
