#!/usr/bin/env python

import tempfile
import os.path
import shutil

import vcr
import numpy as np
from shakemap.coremods.dyfi import _get_dyfi_dataframe
from libcomcat.search import get_event_by_id


def get_datadir():
    # where is this script?
    homedir = os.path.dirname(os.path.abspath(__file__))
    datadir = os.path.join(homedir, '..', '..', 'data')
    return datadir


def test_geocoded():
    # first, test event with 10k and 1k geojson data
    eventid = 'ci14607652'
    datadir = get_datadir()
    tape_file1 = os.path.join(datadir, 'vcr_event1.yaml')

    with vcr.use_cassette(tape_file1):
        detail = get_event_by_id(eventid)
        df, msg = _get_dyfi_dataframe(detail)

    np.testing.assert_almost_equal(df['INTENSITY'].sum(), 4510.1)

    # next, test event with only geocoded (?) resolution text data
    eventid = 'ci14745580'
    tape_file2 = os.path.join(datadir, 'vcr_event2.yaml')

    with vcr.use_cassette(tape_file2):
        detail = get_event_by_id(eventid)
        df, msg = _get_dyfi_dataframe(detail)

    np.testing.assert_almost_equal(df['INTENSITY'].sum(), 800.4)


def test_dyfi():
    eventid = 'nc72282711'
    try:
        tdir = tempfile.mkdtemp()
        datadir = get_datadir()
        tape_file3 = os.path.join(datadir, 'vcr_event3.yaml')

        with vcr.use_cassette(tape_file3):
            detail = get_event_by_id(eventid)
            dataframe, msg = _get_dyfi_dataframe(detail)

    except Exception:
        assert 1 == 2
    finally:
        if os.path.isdir(tdir):
            shutil.rmtree(tdir)

    # Test reading a file
    eventdir = 'eventdata/nc72282711/current/data'
    testfile = os.path.join(datadir,eventdir,'dyfi_geo_10km.geojson')
    dataframe, msg = _get_dyfi_dataframe(None,testfile)
    assert len(dataframe) == 203


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_geocoded()
    test_dyfi()
