#!/usr/bin/env python

# stdlib imports
import os.path
import sys
import io

# third party
import numpy as np
import pytest

# local imports
from shakelib.rupture.origin import Origin


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
month="05" day="12" hour="06" minute="28" second="01" timezone="GMT"
depth="19.0"
locstring="EASTERN SICHUAN, CHINA" created="1211173621" otime="1210573681"
type="" />"""
    source_text = "mech=RS"
    ffile = io.StringIO(fault_text)  # noqa
    efile = io.StringIO(event_text)
    sfile = io.StringIO(source_text)
    origin = Origin.fromFile(efile, sourcefile=sfile)

    testdict = {'mag': 7.9,
                'eventsourcecode': 'us2008ryan',
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

    assert origin.mech == "RS"
    origin.setMechanism("SS")
    assert origin.mech == "SS"
    origin.setMechanism("SS", rake=10)
    assert origin.mech == "SS"
    assert origin.rake == 10.0
    assert origin.dip == 90.0
    origin.setMechanism("SS", rake=-350)
    assert origin.rake == 10.0
    origin.setMechanism("SS", rake=370)
    assert origin.rake == 10.0
    with pytest.raises(Exception) as a:
        origin.setMechanism("SS", dip=100)
    with pytest.raises(Exception) as a:
        origin.setMechanism("Strike slip")
    # Rake too large
    with pytest.raises(Exception) as a:
        origin.setMechanism("SS", rake=1111)
    # Lat is greater than 90
    with pytest.raises(Exception) as a:
        event_text = """<?xml version="1.0" encoding="US-ASCII"
standalone="yes"?> <earthquake id="2008" lat="91.9858" lon="103.3639"
mag="7.9" year="2008" month="05" day="12" hour="06" minute="28" second="01"
timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA"
created="1211173621" otime="1210573681" type="" />"""
        efile = io.StringIO(event_text)
        origin = Origin.fromFile(efile)
    # Lon is greater than 180
    with pytest.raises(Exception) as a:
        event_text = """<?xml version="1.0" encoding="US-ASCII"
standalone="yes"?> <earthquake id="2008" lat="31.9858" lon="183.3639"
mag="7.9" year="2008" month="05" day="12" hour="06" minute="28" second="01"
timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA"
created="1211173621" otime="1210573681" type="" />"""
        efile = io.StringIO(event_text)
        origin = Origin.fromFile(efile)
    # No event id
    with pytest.raises(Exception) as a:
        event_text = """<?xml version="1.0" encoding="US-ASCII"
standalone="yes"?> <earthquake lat="30.9858" lon="103.3639" mag="7.9"
year="2008" month="05" day="12" hour="06" minute="28" second="01"
timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA"
created="1211173621" otime="1210573681" type="" />"""
        efile = io.StringIO(event_text)
        origin = Origin.fromFile(efile)
    # Put mech in event keys
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="2008" lat="30.9858" lon="103.3639" mag="7.9" year="2008"
month="05" day="12" hour="06" minute="28" second="01" timezone="GMT"
depth="19.0" locstring="EASTERN SICHUAN, CHINA" created="1211173621"
otime="1210573681" type="" mech="SS"/>"""
    efile = io.StringIO(event_text)
    origin = Origin.fromFile(efile)
    assert origin.mech == 'SS'
    # Empty mech
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="2008" lat="30.9858" lon="103.3639" mag="7.9" year="2008"
month="05" day="12" hour="06" minute="28" second="01" timezone="GMT"
depth="19.0" locstring="EASTERN SICHUAN, CHINA" created="1211173621"
otime="1210573681" type="" mech=""/>"""
    efile = io.StringIO(event_text)
    origin = Origin.fromFile(efile)
    assert origin.mech == 'ALL'
    # Mech not acceptable value
    with pytest.raises(Exception) as a:  # noqa
        event_text = """<?xml version="1.0" encoding="US-ASCII"
standalone="yes"?> <earthquake id="2008" lat="31.9858" lon="103.3639"
mag="7.9" year="2008" month="05" day="12" hour="06" minute="28" second="01"
timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA"
created="1211173621" otime="1210573681" type="" mech="Strike slip"/>"""
        efile = io.StringIO(event_text)
        origin = Origin.fromFile(efile)


if __name__ == "__main__":
    test_origin()
