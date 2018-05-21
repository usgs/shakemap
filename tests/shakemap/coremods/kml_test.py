#!/usr/bin/env python

import os.path
import tempfile
import shutil
import zipfile

import configobj
from PIL import Image
from lxml import etree

from shakemap.utils.config import get_config_paths
from shakemap.coremods.kml import (create_overlay_image,
                                   create_overlay_kmz,
                                   create_station_kmz,
                                   create_contour_kml,
                                   flip_rgb,
                                   create_styles,
                                   make_placemark)
from shakelib.utils.containers import ShakeMapOutputContainer
from shakemap.utils.utils import path_macro_sub
from impactutils.colors.cpalette import ColorPalette


def test_overlay_img():
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
        im = Image.open(outfile)
        assert (im.format, im.size, im.mode) == ('PNG', (357, 297), 'RGBA')
    except Exception as e:
        assert 1 == 2
    finally:
        shutil.rmtree(tempdir)


def test_overlay_kml():
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
        nlist = myzip.namelist()
        assert sorted(nlist) == ['ii_overlay.png', 'overlay.kml']
        bcomp = b'<kml>\n  <NetworkLinkControl>\n    <minRefreshPeriod>300</minRefreshPeriod>\n  </NetworkLinkControl>\n  <GroundOverlay>\n    <name>Intensity Overlay</name>\n    <color>ffffffff</color>\n    <drawOrder>0</drawOrder>\n    <Icon>\n      <refreshInterval>300</refreshInterval>\n      <refreshMode>onInterval</refreshMode>\n      <href>ii_overlay.png</href>\n    </Icon>\n    <LatLonBox>\n      <north>35.5167</north>\n      <south>33.0500</south>\n      <east>-117.0333</east>\n      <west>-120.0000</west>\n    </LatLonBox>\n  </GroundOverlay>\n</kml>\n'
        mybytes = myzip.read('overlay.kml')
        assert mybytes == bcomp
    except Exception as e:
        assert 1 == 2
    finally:
        os.remove(kmzfile)


def test_station_kml():
    tempdir = tempfile.mkdtemp()
    try:
        homedir = os.path.dirname(os.path.abspath(__file__))
        cfile = os.path.join(homedir, '..', '..', 'data', 'containers',
                             'northridge', 'shake_result.hdf')
        container = ShakeMapOutputContainer.load(cfile)
        install_path, data_path = get_config_paths()
        kmzfile = create_station_kmz(container, tempdir)
        myzip = zipfile.ZipFile(kmzfile, mode='r')

        myzip.close()
    except Exception as e:
        assert 1 == 2
    finally:
        shutil.rmtree(tempdir)


# def test_placemark():
#     # create a dummy station
#     pgv = {'flag': '0',
#            'name': 'pgv',
#            'ln_sigma': 0.0,
#            'units': 'cm/sec',
#            'value': 5.0}
#     pga = {'flag': '0',
#            'name': 'pga',
#            'ln_sigma': 0.0,
#            'units': '%g',
#            'value': 100}
#     psa03 = {'flag': '0',
#              'name': 'sa(0.3)',
#              'ln_sigma': 0.0,
#              'units': '%g',
#              'value': 101}
#     psa10 = {'flag': '0',
#              'name': 'sa(1.0)',
#              'ln_sigma': 0.0,
#              'units': '%g',
#              'value': 102}
#     psa30 = {'flag': '0',
#              'name': 'sa(3.0)',
#              'ln_sigma': 0.0,
#              'units': '%g',
#              'value': 103}
#     amplitudes = [pgv, pga, psa03, psa10, psa30]
#     channels = [{'name': 'HNZ',
#                  'amplitudes': amplitudes}]
#     station = {'id': '1',
#                'geometry': {'type': 'Point', 'coordinates': [-119, 34]},
#                'properties': {'distance': 100,
#                               'intensity': 5.5,
#                               'network': 'US',
#                               'location': 'Somewhere in LA',
#                               'channels': channels}}
#     MY_NAMESPACES = {'gx': 'http://www.google.com/kml/ext/2.2',
#                      'kml': 'http://www.opengis.net/kml/2.2',
#                      'atom': 'http://www.w3.org/2005/Atom',
#                      None:    'http://www.opengis.net/kml/2.2'}
#     root = etree.Element("kml", nsmap=MY_NAMESPACES)
#     nlink = etree.SubElement(root, "NetworkLinkControl")
#     nperiod = etree.SubElement(nlink, "minRefreshPeriod")
#     nperiod.text = '21600'
#     document = etree.SubElement(root, 'Document')
#     create_styles(document)
#     cpalette = ColorPalette.fromPreset('mmi')
#     make_placemark(document, station, cpalette)
#     datadir = os.path.expanduser('~')
#     kmlfile = os.path.join(datadir, 'station.kml')
#     tree = etree.ElementTree(root)
#     tree.write(kmlfile, pretty_print=True)
#     trifile = os.path.join(datadir, 'triangle.png')
#     cirfile = os.path.join(datadir, 'circle.png')
#     kmzfile = os.path.join(datadir, 'station.kmz')
#     station_zip = zipfile.ZipFile(kmzfile,
#                                   mode='w',
#                                   compression=zipfile.ZIP_DEFLATED)
#     station_zip.write(kmlfile, 'station.kml')
#     station_zip.write(trifile, arcname='triangle.png')
#     station_zip.write(cirfile, arcname='circle.png')
#     station_zip.close()


if __name__ == '__main__':
    # test_placemark()
    test_station_kml()
    test_overlay_img()
    test_overlay_kml()
