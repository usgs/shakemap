#!/usr/bin/env python

#stdlib imports
import os.path

#local imports
from shakemap.products.genprod import make_xml_grid

#third party imports
from shakelib.utils.containers import OutputContainer

def test_genprod():
    homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
    datafile = os.path.join(homedir,'..','data','shake_result.hdf')
    container = OutputContainer.load(datafile)
    config = container.getDictionary('config')
    component = config['interp']['component']
    shake_grid = make_xml_grid(container,component)
    assert shake_grid.getLayer('MMI').getData().sum() == 695304.2594600243

if __name__ == '__main__':
    test_genprod()
