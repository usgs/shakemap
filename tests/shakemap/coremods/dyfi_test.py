#!/usr/bin/env python

import tempfile
import os.path
import shutil
import argparse

import numpy as np
from shakemap.coremods.dyfi import _get_dyfi_dataframe
from shakelib.utils.containers import ShakeMapOutputContainer
from libcomcat.search import get_event_by_id

def test_geocoded():
    # first, test event with 10k and 1k geojson data
    eventid = 'ci14607652'
    detail = get_event_by_id(eventid)
    df,msg = _get_dyfi_dataframe(detail)
    np.testing.assert_almost_equal(df['intensity'].sum(),4510.1)
    
    # next, test event with only geocoded (?) resolution text data
    eventid = 'ci14745580'
    detail = get_event_by_id(eventid)
    df,msg = _get_dyfi_dataframe(detail)
    np.testing.assert_almost_equal(df['intensity'].sum(),1162.8)

def test_dyfi():
    eventid = 'nc72282711'
    try:
        tdir = tempfile.mkdtemp()
        detail = get_event_by_id(eventid)
        dataframe,msg = _get_dyfi_dataframe(detail)
    except Exception as e:
        assert 1==2
    finally:
        if os.path.isdir(tdir):
            shutil.rmtree(tdir)
    
if __name__ == '__main__':
    test_geocoded()
    test_dyfi()
