#!/usr/bin/env python

import os.path
import sys

import shakelib.conversions.imt.bommer_alarcon_2006 as ba06

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..', '..'))
sys.path.insert(0, shakedir)


def test_bommeralarcon2006():
    # Inputs
    PGVin = 10
    PSA05in = 0.1

    # BommerAlarcon2006
    PSA05out = ba06.BommerAlarcon2006.pgv2psa05(PGVin)
    PGVout = ba06.BommerAlarcon2006.psa052pgv(PSA05in)
    vfact = ba06.BommerAlarcon2006.getVfact()

    assert abs(PSA05out - 0.2038735983690112) < 0.0001
    assert abs(PGVout - 4.905000) < 0.001
    assert abs(vfact - 49.050000) < 0.001


if __name__ == '__main__':
    test_bommeralarcon2006()
