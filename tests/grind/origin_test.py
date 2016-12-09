#!/usr/bin/env python

# stdlib imports
import os.path
import sys
import io

# third party
import numpy as np
import pytest

# local imports
from shakemap.grind.origin import Origin

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)


def test_origin():
    fault_text = """30.979788       103.454422      1
31.691615       104.419160      1
31.723569       104.374760      1
32.532213       105.220821      1
32.641450       105.135050      20
31.846790       104.246202      20
31.942158       104.205286      20
31.290105       103.284388      20
30.979788       103.454422      1"""
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="2008ryan" lat="30.9858" lon="103.3639" mag="7.9" year="2008" 
month="05" day="12" hour="06" minute="28" second="01" timezone="GMT" depth="19.0" 
locstring="EASTERN SICHUAN, CHINA" created="1211173621" otime="1210573681" type="" />"""
    source_text = "mech=RS"
    ffile = io.StringIO(fault_text)
    efile = io.StringIO(event_text)
    sfile = io.StringIO(source_text)
    origin = Origin.fromFile(efile, sourcefile=sfile)

    testdict = {'mag': 7.9,
                'id': '2008ryan',
                'locstring': 'EASTERN SICHUAN, CHINA',
                'mech': 'RS',
                'lon': 103.3639,
                'lat': 30.9858,
                'depth': 19.0}
    for key in testdict.keys():
        value = eval('origin.%s' % key)
        if type(value) is str:
            assert testdict[key] == value
        if type(value) is float:
            np.testing.assert_almost_equal(testdict[key], value)


if __name__ == "__main__":
    test_origin()
