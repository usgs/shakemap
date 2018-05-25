#!/usr/bin/env python

# Standard library imports
import os.path
import sys

# Third party imports
import pytest

# Local imports
from shakelib.conversions.imt.bommer_alarcon_2006 import BommerAlarcon2006


homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..', '..'))
sys.path.insert(0, shakedir)


def test_bommeralarcon2006():
    # Inputs
    PGVin = 10
    PSA05in = 0.1
    ba06 = BommerAlarcon2006()

    # Test that the correct inputs are returned for valid outputs
    input1 = ba06.getInputIMT('pgv  ')
    input2 = ba06.getInputIMT('PGV')
    input3 = ba06.getInputIMT('psa05  ')
    input4 = ba06.getInputIMT('PSA05')
    assert len(input1) == len(input2) == 1
    assert input1[0] == input2[0] == 'PSA05'
    assert len(input3) == len(input4) == 1
    assert input3[0] == input4[0] == 'PGV'

    # Test that invalid outputs return None
    input5 = ba06.getInputIMT('INVALID')
    assert input5 == None

    # Test valid conversions
    PSA05out = ba06.convertAmps('pgv', 'PSA05', PGVin)
    PGVout = ba06.convertAmps('PSA05', 'PGV', PSA05in)
    vfact = ba06.getConversionFactor()
    assert abs(PSA05out - 0.2038735983690112) < 0.0001
    assert abs(PGVout - 4.905000) < 0.001
    assert abs(vfact - 49.050000) < 0.001

    # Test invalid conversions types
    with pytest.raises(Exception) as a:
        tmp = ba06.convertAmps('INVALID', 'PSA05', PGVin)
    with pytest.raises(Exception) as a:
        tmp = ba06.convertAmps('PGV', 'INVALID', PGVin)
    with pytest.raises(Exception) as a:
        tmp = ba06.convertAmps('INVALID', 'INVALID', PGVin)

if __name__ == '__main__':
    test_bommeralarcon2006()
