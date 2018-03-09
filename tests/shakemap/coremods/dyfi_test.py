#!/usr/bin/env python

import tempfile
import os.path
import shutil
import argparse

import numpy as np
from shakemap.coremods.dyfi import _get_dyfi_dataframe,_dataframe_to_xml
from shakelib.utils.containers import ShakeMapOutputContainer
from libcomcat.search import get_event_by_id

def test_dyfi():
    eventid = 'nc72282711'
    try:
        tdir = tempfile.mkdtemp()
        detail = get_event_by_id(eventid)
        dataframe,msg = _get_dyfi_dataframe(detail)
        _dataframe_to_xml(dataframe,eventid,tdir,reference='DYFI System')
    except Exception as e:
        assert 1==2
    finally:
        if os.path.isdir(tdir):
            shutil.rmtree(tdir)
    
if __name__ == '__main__':
    test_dyfi()
