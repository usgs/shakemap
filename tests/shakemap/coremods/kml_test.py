#!/usr/bin/env python

import os.path
import tempfile
import shutil
import zipfile
import logging
from xml.dom import minidom

from shakemap.coremods.base import Contents
from shakemap.utils.config import (get_config_paths)
from shakemap.coremods.kml import (create_kmz)
from impactutils.io.smcontainers import ShakeMapOutputContainer


def test_create_kmz():
    tempdir = tempfile.mkdtemp()
    try:
        homedir = os.path.dirname(os.path.abspath(__file__))
        cfile = os.path.join(homedir, '..', '..', 'data', 'containers',
                             'northridge', 'shake_result.hdf')
        container = ShakeMapOutputContainer.load(cfile)
        install_path, data_path = get_config_paths()

        logger = logging.getLogger(__name__)
        contents = Contents(None, None, 'northridge')
        kmzfile = create_kmz(container, tempdir, logger, contents)
        myzip = zipfile.ZipFile(kmzfile, mode='r')
        kmlstr = myzip.read('shakemap.kml').decode('utf-8')
        root = minidom.parseString(kmlstr)
        document = root.getElementsByTagName('Document')[0]
        folders = document.getElementsByTagName('Folder')
        names = []
        nstations = 0
        nmmi = 0
        for folder in folders:
            name = folder.getElementsByTagName('name')[0].firstChild.data
            names.append(name)
            if name == 'Instrumented Stations':
                nstations = len(folder.getElementsByTagName('Placemark'))
            elif name == 'Macroseismic Stations':
                nmmi = len(folder.getElementsByTagName('Placemark'))
        assert sorted(names) == ['Contours', 'Instrumented Stations',
                                 'MMI 4 Polygons', 'MMI 5 Polygons',
                                 'MMI 6 Polygons', 'MMI 7 Polygons',
                                 'MMI 8 Polygons', 'MMI 8.5 Polygons',
                                 'MMI Contours', 'MMI Labels',
                                 'MMI Polygons', 'Macroseismic Stations',
                                 'PGA Contours', 'PGV Contours',
                                 'SA(0.3) Contours', 'SA(1.0) Contours',
                                 'SA(3.0) Contours']
        assert nstations == 185
        assert nmmi == 547
        myzip.close()

    except Exception as e:
        print(str(e))
        assert 1 == 2
    finally:
        shutil.rmtree(tempdir)


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_create_kmz()
