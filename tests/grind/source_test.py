#!/usr/bin/env python

# stdlib imports
import os.path
import sys
import io

# hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
# put this at the front of the system path, ignoring any installed
# shakemap stuff
sys.path.insert(0, shakedir)

# third party
from openquake.hazardlib.gsim import abrahamson_2014
import numpy as np
import pytest

# local imports
from shakemap.grind.source import Source


def test_source():
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
<earthquake id="2008ryan" lat="30.9858" lon="103.3639" mag="7.9" year="2008" month="05" day="12" hour="06" minute="28" second="01" timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA" created="1211173621" otime="1210573681" type="" />
    """
    source_text = "mech=RS"
    ffile = io.StringIO(fault_text)
    efile = io.StringIO(event_text)
    sfile = io.StringIO(source_text)
    source = Source.readFromFile(efile, faultfile=ffile, sourcefile=sfile)

    gmpe = abrahamson_2014.AbrahamsonEtAl2014()
    rupture = source.getRuptureContext([gmpe])
    testdict = {'mag': 7.9,
                'strike': -133.083550974,
                'dip': 49.8524115024,
                'rake': 45.0,
                'ztor': 0.999999999995,
                'hypo_lon': 103.3639,
                'hypo_lat': 30.9858,
                'hypo_depth': 19.0,
                'width': 27.8623813381}
    for key in testdict.keys():
        value = eval('rupture.%s' % key)
        np.testing.assert_almost_equal(testdict[key], value)

    mech = 'RS'
    exp_dip = 40
    exp_rake = 90
    source.setMechanism(mech)
    assert source.getEventParam('dip') == exp_dip
    assert source.getEventParam('rake') == exp_rake
    source.setMechanism('ALL', dip=45, rake=315)
    assert source.getEventParam('rake') == -45

    # this should raise an exception
    with pytest.raises(Exception) as e_info:
        source.setMechanism('ALL', dip=110)
    with pytest.raises(Exception) as e_info:
        source.setMechanism('ALL', rake=620)
