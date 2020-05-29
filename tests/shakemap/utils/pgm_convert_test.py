#!/usr/bin/env python
import pathlib
from shakemap.utils.pgm_convert import get_stations_dict
from numpy.testing import assert_almost_equal


def test_get_stations_dict():
    thisdir = pathlib.Path(__file__).parent.absolute()
    datadir = thisdir.parents[1] / 'data' / 'ampdata'

    mmifile = datadir / 'ceresis_amps.xlsx'
    mmidict, mmirows = get_stations_dict(mmifile)
    assert len(mmidict['features']) == 21
    assert mmidict['features'][0]['properties']['intensity'] == 8

    instfile = datadir / 'mem_ridgecrest.xlsx'
    instdict, instrows = get_stations_dict(instfile)
    assert len(instdict['features']) == 358
    channel1 = instdict['features'][0]['properties']['channels'][0]['name']
    assert channel1 == 'H1'
    value1 = (instdict['features'][0]['properties']['channels']
              [0]['amplitudes'][0]['value'])

    assert_almost_equal(value1, 2.243548)


if __name__ == '__main__':
    test_get_stations_dict()
