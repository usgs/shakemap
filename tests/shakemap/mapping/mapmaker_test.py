#!/usr/bin/env python

import os.path

from shakelib.utils.containers import ShakeMapOutputContainer
from shakemap.mapping.mapmaker import draw_intensity


def old_test_mapmaker():
    homedir = os.path.dirname(os.path.abspath(
        __file__))  # where is this script?
    shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))
    out_file = os.path.join(shakedir, 'tests', 'data',
                            'containers', 'northridge',
                            'shake_result.hdf')
    container = ShakeMapOutputContainer.load(out_file)
    topofile = os.path.join(homedir, '..', '..', 'install', 'data',
                            'mapping', 'CA_topo.grd')
    png, cities = draw_intensity(container, topofile,
                                 oceanfile, oceangridfile,
                                 cityfile, basename)


if __name__ == '__main__':
    test_mapmaker()
