#!/usr/bin/env python

import tempfile
import os.path

import numpy as np
from shakemap.coremods.transfer import _transfer
from shakelib.utils.containers import ShakeMapOutputContainer

def test_transfer():
    homedir = os.path.dirname(os.path.abspath(__file__))
    cfile = os.path.join(homedir,'..','..','data/containers/northridge/shake_result.hdf')
    products_dir = os.path.join(homedir,'..','..','/data/eventdata/northridge/current/')
    container = ShakeMapOutputContainer.load(cfile)
    try:
        tdir = tempfile.mkdtemp()
        config = {'copy':{'local':{'remote_directory':tdir}}}
        _transfer(config,container,products_dir)
    except Exception as e:
        pass
    finally:
        if os.path.isdir(tdir):
            shutil.rmtree(tdir)

if __name__ == '__main__':
    test_transfer()
