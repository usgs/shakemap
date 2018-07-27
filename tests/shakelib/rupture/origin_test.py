#!/usr/bin/env python

# stdlib imports
import io
import os
import tempfile

# third party
import numpy as np
import pytest

# local imports
from shakelib.rupture.origin import Origin
from shakelib.rupture.origin import write_event_file
from impactutils.time.ancient_time import HistoricTime
from shakelib.rupture import constants


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
<earthquake id="2008ryan" lat="30.9858" lon="103.3639" mag="7.9"
time="2008-05-12T06:28:01Z"
depth="19.0" netid="us" network=""
locstring="EASTERN SICHUAN, CHINA"
mech="" />"""
    source_text = "mech=RS"
    ffile = io.StringIO(fault_text)  # noqa
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
mag="7.9" time="2008-05-12T06:28:01Z"
timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA"
created="1211173621" otime="1210573681" mech="" netid="us" network=""/>"""
        efile = io.StringIO(event_text)
        origin = Origin.fromFile(efile)
    # Lon is greater than 180
    with pytest.raises(Exception) as a:
        event_text = """<?xml version="1.0" encoding="US-ASCII"
standalone="yes"?> <earthquake id="2008" lat="31.9858" lon="183.3639"
mag="7.9" time="2008-05-12T06:28:01Z"
timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA"
created="1211173621" otime="1210573681" type="" netid="us" network=""/>"""
        efile = io.StringIO(event_text)
        origin = Origin.fromFile(efile)
    # No event id
    with pytest.raises(Exception) as a:
        event_text = """<?xml version="1.0" encoding="US-ASCII"
standalone="yes"?> <earthquake lat="30.9858" lon="103.3639" mag="7.9"
time="2008-05-12T06:28:01Z"
timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA"
created="1211173621" otime="1210573681" type="" netid="us" network=""/>"""
        efile = io.StringIO(event_text)
        origin = Origin.fromFile(efile)
    # Put mech in event keys
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="2008" lat="30.9858" lon="103.3639" mag="7.9"
time="2008-05-12T06:28:01Z"
depth="19.0" locstring="EASTERN SICHUAN, CHINA" created="1211173621"
otime="1210573681" type="" mech="SS" netid="us" network=""/>"""
    efile = io.StringIO(event_text)
    origin = Origin.fromFile(efile)
    assert origin.mech == 'SS'
    # Empty mech
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="2008" lat="30.9858" lon="103.3639" mag="7.9"
time="2008-05-12T06:28:01Z"
depth="19.0" locstring="EASTERN SICHUAN, CHINA" created="1211173621"
otime="1210573681" type="" mech="" netid="us" network=""/>"""
    efile = io.StringIO(event_text)
    origin = Origin.fromFile(efile)
    assert origin.mech == 'ALL'
    # Mech not acceptable value
    with pytest.raises(Exception) as a:  # noqa
        event_text = """<?xml version="1.0" encoding="US-ASCII"
standalone="yes"?> <earthquake id="2008" lat="31.9858" lon="103.3639"
mag="7.9"
time="2008-05-12T06:28:01Z"
depth="19.0" locstring="EASTERN SICHUAN, CHINA"
created="1211173621" otime="1210573681" type="" mech="Strike slip"
netid="us" network=""/>"""
        efile = io.StringIO(event_text)
        origin = Origin.fromFile(efile)

    # Missing keys
    with pytest.raises(KeyError):
        event_text = """<?xml version="1.0" encoding="US-ASCII"
standalone="yes"?> <earthquake id="2008"
mag="7.9"
time="2008-05-12T06:28:01Z"
depth="19.0" locstring="EASTERN SICHUAN, CHINA"
created="1211173621" otime="1210573681" type=""
network=""/>"""
        efile = io.StringIO(event_text)
        origin = Origin.fromFile(efile)

    # Use "type" instead of "mech"
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="2008ryan" lat="30.9858" lon="103.3639" mag="7.9"
time="2008-05-12T06:28:01Z"
depth="19.0" netid="us" network=""
locstring="EASTERN SICHUAN, CHINA"
type="RS" />"""
    efile = io.StringIO(event_text)
    origin = Origin.fromFile(efile)
    assert origin.mech == 'RS'

    # No rake or mech
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="2008ryan" lat="30.9858" lon="103.3639" mag="7.9"
time="2008-05-12T06:28:01Z"
depth="19.0" netid="us" network=""
locstring="EASTERN SICHUAN, CHINA"
reference="Smith, et al. (2019)"
 />"""
    efile = io.StringIO(event_text)
    origin = Origin.fromFile(efile)
    assert origin.rake == 0.0

    hypo = origin.getHypo()
    assert hypo.x == origin.lon
    assert hypo.y == origin.lat
    assert hypo.z == origin.depth

    # Write the origin to a file
    event = {}
    event['id'] = 'us2000ryan'
    event['netid'] = 'us'
    event['network'] = 'USGS Network'
    event['lat'] = 30.9858
    event['lon'] = 103.3639
    event['depth'] = 19.0
    event['mag'] = 7.9
    event['time'] = HistoricTime.strptime('2008-05-12T06:28:01.0Z',
                                          constants.TIMEFMT)
    event['locstring'] = "EASTERN SICHUAN, CHINA"
    event['mech'] = 'RS'
    event['reference'] = "Smith, et al. (2019)"
    event['productcode'] = "us2000ryan"

    tfile = tempfile.NamedTemporaryFile()
    xmlfile = tfile.name
    tfile.close()
    res = write_event_file(event, xmlfile)
    if res is not None:
        print(res)
        assert False

    origin = Origin.fromFile(xmlfile)
    os.remove(xmlfile)
    assert origin.id == event['id']
    assert origin.netid == event['netid']
    assert origin.network == event['network']
    assert origin.time == event['time']


if __name__ == "__main__":
    test_origin()
