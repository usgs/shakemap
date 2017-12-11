#!/usr/bin/env python


import numpy as np

import os.path
import sys

from shakelib.conversions.imt.newmark_hall_1982 import NewmarkHall1982

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..', '..'))
sys.path.insert(0, shakedir)


def test_newmarkhall1982():
    # Inputs
    # PGVin = np.log(10)
    PSA10in = np.log(0.1)
    sd = 0.6

    # NewmarkHall1982
    PGVout, PGVsdout = NewmarkHall1982.psa102pgv(PSA10in, sd)
    mfact = NewmarkHall1982.getConversionFactor()
    lnsig = NewmarkHall1982.getLnSigma()

    assert abs(PGVout - np.log(9.46658)) < 0.001
    assert abs(PGVsdout - 0.790489) < 0.001
    assert abs(mfact - 94.6658) < 0.001
    assert abs(lnsig - 0.5146578) < 0.001


if __name__ == '__main__':
    test_newmarkhall1982()
