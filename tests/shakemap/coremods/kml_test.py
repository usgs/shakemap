#!/usr/bin/env python

import os.path
import tempfile
import shutil
import zipfile

import configobj

from shakemap.utils.config import get_config_paths
from shakemap.coremods.kml import create_overlay_image, create_overlay_kml
from shakelib.utils.containers import ShakeMapOutputContainer
from shakemap.utils.utils import path_macro_sub


def test_overlay():
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
        outfile = os.path.join(tempdir, 'test.png')
        create_overlay_image(container, oceanfile, outfile)
    except Exception as e:
        assert 1 == 2
    finally:
        shutil.rmtree(tempdir)

    # TODO - how to test image is valid?


def test_kml():
    homedir = os.path.dirname(os.path.abspath(__file__))
    cfile = os.path.join(homedir, '..', '..', 'data', 'containers',
                         'northridge', 'shake_result.hdf')
    container = ShakeMapOutputContainer.load(cfile)
    install_path, data_path = get_config_paths()
    # datadir = os.path.join(data_path, eventid, 'current', 'products')
    datadir = os.path.join(os.path.expanduser('~'))
    product_config_file = os.path.join(
        install_path, 'config', 'products.conf')
    pconfig = configobj.ConfigObj(product_config_file)
    oceanfile = pconfig['products']['mapping']['layers']['lowres_oceans']
    oceanfile = path_macro_sub(oceanfile, install_path, data_path)
    try:
        kmzfile = create_overlay_kml(container, oceanfile, datadir)
        # check the kmz file contents
        myzip = zipfile.ZipFile(kmzfile)

    except Exception as e:
        assert 1 == 2
    finally:
        os.remove(kmzfile)


if __name__ == '__main__':
    test_overlay()
    test_kml()
