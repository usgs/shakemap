#!/usr/bin/env python

import tempfile
import os.path
import shutil

import numpy as np
from shakemap.coremods.transfer import TransferModule
from shakelib.utils.containers import ShakeMapOutputContainer

def test_transfer():
    homedir = os.path.dirname(os.path.abspath(__file__))
    cfile = os.path.join(homedir,'..','..','data','containers',
                         'northridge','shake_result.hdf')
    products_dir = os.path.join(homedir,'..','..','data',
                                'eventdata','northridge','current')
    eventid = 'ci3144585'
    
    container = ShakeMapOutputContainer.load(cfile)
    try:
        tdir = tempfile.mkdtemp()
        remote_dir = os.path.join(tdir,eventid)
        config = {'copy':{'local':{'remote_directory':tdir}}}
        transfermod = TransferModule(eventid)
        transfermod._transfer(config,container,products_dir)
        nfiles = len(os.listdir(remote_dir))
        nsrcfiles = len(os.listdir(products_dir))
        assert nfiles == nsrcfiles
    except Exception as e:
        pass
    finally:
        if os.path.isdir(tdir):
            shutil.rmtree(tdir)

if __name__ == '__main__':
    test_transfer()
