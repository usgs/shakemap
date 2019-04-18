#!/usr/bin/env python

import os.path
from tempfile import mkdtemp
import shutil
import copy

from impactutils.io.smcontainers import ShakeMapOutputContainer
from impactutils.mapping.city import Cities
from shakemap.mapping.mapmaker import draw_map
from shakemap.coremods.mapping import get_text_strings
from shakemap.utils.config import get_data_path

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

    outpath = mkdtemp()

    model_config = container.getConfig()
    comp = container.getComponents('MMI')[0]
    textfile = os.path.join(get_data_path(), 'mapping',
                            'map_strings.en')
    text_dict = get_text_strings(textfile)

    cities = Cities.fromDefault()
    d = {'imtype': 'MMI',
         'topogrid': topogrid,
         'allcities': cities,
         'states_provinces': None,
         'countries': None,
         'oceans': None,
         'lakes': None,
         'roads': None,
         'faults': None,
         'datadir': outpath,
         'operator': 'NEIC',
         'filter_size': 10,
         'info': info,
         'component': comp,
         'imtdict': container.getIMTGrids('MMI', comp),
         'ruptdict': copy.deepcopy(container.getRuptureDict()),
         'stationdict': container.getStationDict(),
         'config': model_config,
         'tdict': text_dict
         }

    try:
        fig1, fig2 = draw_map(d)
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

    outpath = mkdtemp()
    filter_size = 10
    model_config = container.getConfig()
    comp = container.getComponents('PGA')[0]
    textfile = os.path.join(get_data_path(), 'mapping',
                            'map_strings.en')
    text_dict = get_text_strings(textfile)

    cities = Cities.fromDefault()
    d = {'imtype': 'PGA',
         'topogrid': topogrid,
         'allcities': cities,
         'states_provinces': None,
         'countries': None,
         'oceans': None,
         'lakes': None,
         'roads': None,
         'faults': None,
         'datadir': outpath,
         'operator': 'NEIC',
         'filter_size': filter_size,
         'info': info,
         'component': comp,
         'imtdict': container.getIMTGrids('PGA', comp),
         'ruptdict': copy.deepcopy(container.getRuptureDict()),
         'stationdict': container.getStationDict(),
         'config': model_config,
         'tdict': text_dict
         }
    try:
        _ = draw_map(d)
    except Exception:
        assert 1 == 2
    finally:
        shutil.rmtree(outpath)


if __name__ == '__main__':
    test_mapmaker_contour()
    test_mapmaker_intensity()
