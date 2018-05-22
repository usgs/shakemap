#!/usr/bin/env python

import os.path
import tempfile
import shutil
import zipfile
import logging
from xml.dom import minidom

import configobj
from PIL import Image
from lxml import etree

from shakemap.utils.config import get_config_paths
from shakemap.coremods.kml import (create_kmz)
from shakelib.utils.containers import ShakeMapOutputContainer
from shakemap.utils.utils import path_macro_sub
from impactutils.colors.cpalette import ColorPalette


def test_create_kmz():
    tempdir = tempfile.mkdtemp()
    try:
        homedir = os.path.dirname(os.path.abspath(__file__))
        cfile = os.path.join(homedir, '..', '..', 'data', 'containers',
                             'northridge', 'shake_result.hdf')
        container = ShakeMapOutputContainer.load(cfile)
        install_path, data_path = get_config_paths()
        product_config_file = os.path.join(
            install_path, 'config', 'products.conf')
        pconfig = configobj.ConfigObj(product_config_file)
        oceanfile = pconfig['products']['mapping']['layers']['lowres_oceans']
        oceanfile = path_macro_sub(oceanfile, install_path, data_path)

        logger = logging.getLogger()

        kmzfile = create_kmz(container, tempdir, oceanfile, logger)
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
        assert sorted(names) == ['Contours',
                                 'Instrumented Stations',
                                 'Macroseismic Stations']
        assert nstations == 185
        assert nmmi == 977
        myzip.close()

    except Exception as e:
        assert 1 == 2
    finally:
        shutil.rmtree(tempdir)


if __name__ == '__main__':
    test_create_kmz()
