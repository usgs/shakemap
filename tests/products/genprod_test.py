#!/usr/bin/env python

#stdlib imports
import os.path

#local imports
from shakemap.products.genprod import make_shake_grid

#third party imports
from shakelib.utils.container import OutputContainer

def test_genprod():
    homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
    datafile = os.path.join(homedir,'..','data','shake_result.hdf')
    container = OutputContainer.loadFromHDF(datafile)
    shake_grid = make_shake_grid(container)

if __name__ == '__main__':
    test_genprod()
