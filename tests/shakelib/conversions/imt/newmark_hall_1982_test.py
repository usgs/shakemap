#!/usr/bin/env python

# Standard library imports
import os.path
import sys

# Third party imports
import numpy as np
import pytest

# Local imports
from shakelib.conversions.imt.newmark_hall_1982 import NewmarkHall1982


homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..', '..'))
sys.path.insert(0, shakedir)


def test_newmarkhall1982():
    # Inputs
    PSA10in = np.log(0.1)
    sd = 0.6
    nh82 = NewmarkHall1982()

    # Test that the correct inputs are returned for valid outputs
    input1 = nh82.getInputIMT('pgv  ')
    input2 = nh82.getInputIMT('PGV')
    assert len(input1) == len(input2) == 1
    assert input1[0] == input2[0] == 'PSA10'

    # Test that invalid outputs return None
    input3 = nh82.getInputIMT('INVALID')
    assert input3 == None

    # Test valid conversions
    PGVout = nh82.convertAmps('psa10', 'PGV', PSA10in)
    PGVsdout = nh82.convertSigmas('psa10', 'PGV', sd)
    mfact = nh82.getConversionFactor()
    lnsig = nh82.getLnSigma()
    assert abs(PGVout - np.log(9.46658)) < 0.001
    assert abs(PGVsdout - 0.790489) < 0.001
    assert abs(mfact - 94.6658) < 0.001
    assert abs(lnsig - 0.5146578) < 0.001

    # Test invalid conversions
    with pytest.raises(Exception) as a:
        tmp = nh82.convertAmps('INVALID', 'PSA05', PGVin)
    with pytest.raises(Exception) as a:
        tmp = nh82.convertAmps('PGV', 'INVALID', PGVin)
    with pytest.raises(Exception) as a:
        tmp = nh82.convertAmps('INVALID', 'INVALID', PGVin)


if __name__ == '__main__':
    test_newmarkhall1982()
