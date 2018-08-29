#!/usr/bin/env python

import os.path
from tempfile import mkdtemp
import shutil

from impactutils.io.smcontainers import ShakeMapOutputContainer
from shakemap.mapping.mapmaker import draw_intensity, draw_contour


def test_mapmaker_intensity():
    homedir = os.path.dirname(os.path.abspath(
        __file__))  # where is this script?
    shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))
    out_file = os.path.join(shakedir, 'tests', 'data',
                            'containers', 'northridge',
                            'shake_result.hdf')
    container = ShakeMapOutputContainer.load(out_file)
    topofile = os.path.join(homedir, '..', '..', 'data', 'install', 'data',
                            'mapping', 'CA_topo.grd')
    oceanfile = os.path.join(homedir, '..', '..', 'data', 'install', 'data',
                             'mapping', 'northridge_ocean.json')
    tdir = mkdtemp()
    basename = os.path.join(tdir, 'testmap')
    basename = os.path.join(os.path.expanduser('~'), 'testmap')
    try:
        pdf, png, cities = draw_intensity(container, topofile,
                                          oceanfile, basename, 'NEIC')
        print(pdf)
    except Exception:
        assert 1 == 2
    finally:
        shutil.rmtree(tdir)


def test_mapmaker_contour():
    homedir = os.path.dirname(os.path.abspath(
        __file__))  # where is this script?
    shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))
    out_file = os.path.join(shakedir, 'tests', 'data',
                            'containers', 'northridge',
                            'shake_result.hdf')
    container = ShakeMapOutputContainer.load(out_file)
    topofile = os.path.join(homedir, '..', '..', 'data', 'install', 'data',
                            'mapping', 'CA_topo.grd')
    oceanfile = os.path.join(homedir, '..', '..', 'data', 'install', 'data',
                             'mapping', 'northridge_ocean.json')
    tdir = mkdtemp()
    # basename = os.path.join(tdir, 'testmap')
    basename = os.path.join(os.path.expanduser('~'), 'testmap')
    filter_size = 10
    try:
        pdf = draw_contour(container, 'PGA', topofile, oceanfile, basename,
                           'NEIC', filter_size, borderfile=None,
                           is_scenario=False)
        print(pdf)
        print()
    except Exception as e:
        print(e)
        assert 1 == 2
    finally:
        shutil.rmtree(tdir)


if __name__ == '__main__':
    test_mapmaker_contour()
    test_mapmaker_intensity()
