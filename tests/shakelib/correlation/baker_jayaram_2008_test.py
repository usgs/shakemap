
import numpy as np

from shakelib.correlation.baker_jayaram_2008 import baker_jayaram_correlation


def test_baker_jayaram_2008():
    pp = np.logspace(-2.0, 1.0, 100)

    outgrid = np.full((100, 100), np.nan)
    for i, pp1 in enumerate(pp):
        for j, pp2 in enumerate(pp):
            outgrid[i, j] = baker_jayaram_correlation(pp1, pp2)

    #
    # This is some crap testing right here
    #
    assert(np.max(outgrid) == 1)
    assert(np.min(outgrid) >= 0)
    assert((outgrid.transpose() == outgrid).all())
