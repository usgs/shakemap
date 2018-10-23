#!/usr/bin/env python

import os.path
from tempfile import mkdtemp
import shutil

from impactutils.io.smcontainers import ShakeMapOutputContainer
from shakemap.mapping.mapmaker import draw_intensity, draw_contour

from mapio.gmt import GMTGrid
from mapio.geodict import GeoDict


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

    info = container.getMetadata()
    xmin = info['output']['map_information']['min']['longitude']
    xmax = info['output']['map_information']['max']['longitude']
    ymin = info['output']['map_information']['min']['latitude']
    ymax = info['output']['map_information']['max']['latitude']
    xmin = float(xmin) - 0.1
    xmax = float(xmax) + 0.1
    ymin = float(ymin) - 0.1
    ymax = float(ymax) + 0.1
    dy = float(info['output']['map_information']
               ['grid_spacing']['latitude'])
    dx = float(info['output']['map_information']
               ['grid_spacing']['longitude'])
    sampledict = GeoDict.createDictFromBox(xmin, xmax, ymin, ymax, dx, dy)
    topogrid = GMTGrid.load(topofile,
                            samplegeodict=sampledict,
                            resample=False)

    oceanfile = os.path.join(homedir, '..', '..', 'data', 'install', 'data',
                             'mapping', 'northridge_ocean.json')
    outpath = mkdtemp()

    try:
        pdf, png, legend = draw_intensity(container, topogrid, oceanfile,
                                          outpath, 'NEIC')
        print(pdf)
    except Exception:
        assert 1 == 2
    finally:
        shutil.rmtree(outpath)


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

    info = container.getMetadata()
    xmin = info['output']['map_information']['min']['longitude']
    xmax = info['output']['map_information']['max']['longitude']
    ymin = info['output']['map_information']['min']['latitude']
    ymax = info['output']['map_information']['max']['latitude']
    xmin = float(xmin) - 0.1
    xmax = float(xmax) + 0.1
    ymin = float(ymin) - 0.1
    ymax = float(ymax) + 0.1
    dy = float(info['output']['map_information']
               ['grid_spacing']['latitude'])
    dx = float(info['output']['map_information']
               ['grid_spacing']['longitude'])
    sampledict = GeoDict.createDictFromBox(xmin, xmax, ymin, ymax, dx, dy)
    topogrid = GMTGrid.load(topofile,
                            samplegeodict=sampledict,
                            resample=False)

    oceanfile = os.path.join(homedir, '..', '..', 'data', 'install', 'data',
                             'mapping', 'northridge_ocean.json')
    outpath = mkdtemp()
    filter_size = 10
    try:
        pdf, png = draw_contour(container, 'PGA', topogrid, oceanfile,
                                outpath, 'NEIC', filter_size)
        print(pdf)
    except Exception:
        assert 1 == 2
    finally:
        shutil.rmtree(outpath)


if __name__ == '__main__':
    test_mapmaker_contour()
    test_mapmaker_intensity()
