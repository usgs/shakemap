#!/usr/bin/env python

import os.path

from configobj import ConfigObj
from validate import Validator
from shakemap.coremods.assemble import (_get_grids)
import numpy as np


def test_assemble_sim():
    homedir = os.path.dirname(os.path.abspath(__file__))
    cfgfile = os.path.join(homedir, '..', '..', 'data', 'simulation',
                           'simulation.conf')
    specfile = os.path.join(homedir, '..', '..', '..', 'shakemap',
                            'data', 'simulationspec.conf')
    simfile = os.path.join(homedir, '..', '..', 'data', 'simulation',
                           'planet9.csv')
    cfgfile = os.path.abspath(cfgfile)
    specfile = os.path.abspath(specfile)
    simfile = os.path.abspath(simfile)
    config = ConfigObj(cfgfile, configspec=specfile)
    vtor = Validator()
    results = config.validate(vtor)
    assert results
    imtgrids = _get_grids(config, simfile)
    pgasum = np.nansum(imtgrids['PGA'].getData())
    pgacomp = -2195.342762128482
    np.testing.assert_almost_equal(pgasum, pgacomp)
    imtlist = list(imtgrids.keys())
    imtcmp = ['PGA', 'PGV', 'SA(0.3)', 'SA(1.0)', 'SA(3.0)']
    assert sorted(imtlist) == sorted(imtcmp)


if __name__ == '__main__':
    test_assemble_sim()
