#!/usr/bin/env python

import os
import os.path
import tempfile
import zipfile
import json
import glob
import subprocess

import numpy as np
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.model import ModelModule
from shakemap.coremods.assemble import AssembleModule
from shakemap.coremods.history import HistoryModule
from shakemap.coremods.contour import ContourModule
from shakemap.coremods.info import InfoModule
from shakemap.coremods.raster import RasterModule
from shakemap.coremods.rupture import RuptureModule
from shakemap.coremods.stations import StationModule
from shakemap.coremods.mapping import MappingModule
from shakemap.coremods.plotregr import PlotRegr
from shakemap.coremods.kml import KMLModule
from shakemap.coremods.gridxml import GridXMLModule, _oq_to_gridxml
from shakemap.coremods.shape import ShapeModule
from impactutils.io.smcontainers import ShakeMapOutputContainer
from shakelib.utils.imt_string import oq_to_file
from mapio.shake import ShakeGrid


#
# Returns a functions that rounds all the elements of a numpy array to a
# certain number of significant digits
# Using a closure here because pep8 complains about assigning a lambda
#
def rounder(digits):
    fmt = '%%.%dg' % digits

    def doround(x):
        return float(fmt % x)

    return np.vectorize(doround)


#
# Order a json object for comparison with another
#
def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


#
# This checks the common failure modes of coremods: the event isn't
# there or the shake_results file hasn't been generated yet
#
def check_failures(evid, datapath, module):

    # Non-existent event (should fail)
    mod = module('not_an_event')
    with pytest.raises(NotADirectoryError):
        mod.execute()

    # Would succeed but we remove event.xml (should fail)
    event_file = os.path.join(datapath, evid, 'current', 'products',
                              'shake_result.hdf')
    os.rename(event_file, event_file + '_safe')
    try:
        mod = module(evid)
        with pytest.raises(FileNotFoundError):
            mod.execute()
    finally:
        os.rename(event_file + '_safe', event_file)


#
# Check the rupture JSON against the input fault
# TODO: Could check more aspects of the JSON
#
def do_rupture(evid, datapath, oc):

    check_failures(evid, datapath, RuptureModule)

    mod = RuptureModule(evid)
    mod.execute()
    mod.writeContents()

    ifile = os.path.join(datapath, evid, 'current', 'products', 'rupture.json')
    with open(ifile, 'r') as infile:
        rupj = json.load(infile)
    ffile = os.path.join(datapath, evid, 'current', 'boat_fault.txt')
    fcoords = []
    with open(ffile, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            lat, lon, depth = line.split()
            fcoords.append((float(lat), float(lon), float(depth)))
    jcoords = rupj['features'][0]['geometry']['coordinates'][0][0]
    for ix, ftup in enumerate(fcoords):
        assert np.allclose(ftup, jcoords[ix])


#
# Checks the info.json file against what is in the output container.
# Not sure how these could be different, but there you go.
#
def do_info(evid, datapath, oc):

    check_failures(evid, datapath, InfoModule)

    mod = InfoModule(evid)
    mod.execute()
    mod.writeContents()

    ifile = os.path.join(datapath, evid, 'current', 'products', 'info.json')
    with open(ifile, 'r') as infile:
        fjson = json.loads(infile.read())
        cjson = oc.getMetadata()
        assert fjson == cjson


#
# Checks the raster grids against the grids in the output container
# TODO: Should also probably check coordinates, resolution, etc.
#
def do_raster(evid, datapath, oc):

    check_failures(evid, datapath, RasterModule)

    mod = RasterModule(evid)
    mod.execute()
    mod.writeContents()

    driver = gdal.GetDriverByName('ENVI')
    driver.Register()

    imts = oc.getIMTs()

    rzip = os.path.join(datapath, evid, 'current', 'products', 'raster.zip')
    with tempfile.TemporaryDirectory() as tmpdirname:
        with zipfile.ZipFile(rzip, 'r') as zip_ref:
            zip_ref.extractall(tmpdirname)
            for imt in imts:
                component, imt = imt.split('/')
                fname = oq_to_file(imt)
                fname = os.path.join(tmpdirname, fname + '_mean.flt')
                rin = gdal.Open(fname, GA_ReadOnly)
                if rin is None:
                    raise RuntimeError("Couldn't open %s" % fname)
                cols = rin.RasterXSize
                rows = rin.RasterYSize
                band = rin.GetRasterBand(1)
                rgrid = band.ReadAsArray(0, 0, cols, rows)
                comp = oc.getComponents(imt)
                cdata = oc.getIMTGrids(imt, comp[0])['mean']
                assert np.allclose(cdata, rgrid)

                fname = oq_to_file(imt)
                fname = os.path.join(tmpdirname, fname + '_std.flt')
                rin = gdal.Open(fname, GA_ReadOnly)
                if rin is None:
                    raise RuntimeError("Couldn't open %s" % fname)
                cols = rin.RasterXSize
                rows = rin.RasterYSize
                band = rin.GetRasterBand(1)
                rgrid = band.ReadAsArray(0, 0, cols, rows)
                comp = oc.getComponents(imt)
                cdata = oc.getIMTGrids(imt, comp[0])['std']
                assert np.allclose(cdata, rgrid)


#
# Checks the grid.xml and uncertainty.xml against the grids in the
# output container.
#
def do_gridxml(evid, datapath, oc):

    check_failures(evid, datapath, GridXMLModule)

    mod = GridXMLModule(evid)
    mod.execute()
    mod.writeContents()

    #
    # Test that the grid.xml grids actually match what's in
    # shake_results.hdf
    #
    imts = oc.getIMTs()

    gxml = os.path.join(datapath, evid, 'current', 'products', 'grid.xml')
    g2d = ShakeGrid.load(gxml)
    layers = g2d.getData()
    for imt in imts:
        component, imt = imt.split('/')
        comp = oc.getComponents(imt)
        cdata = oc.getIMTGrids(imt, comp[0])['mean']
        #
        # Do the same conversion to the container data as is
        # done to the file data
        #
        digits = oc.getIMTGrids(imt, comp[0])['mean_metadata']['digits']
        vfunc = rounder(digits)
        if imt == 'MMI':
            cdata = vfunc(cdata)
        elif imt == 'PGV':
            cdata = vfunc(np.exp(cdata))
        else:
            cdata = vfunc(100 * np.exp(cdata))
        lname = _oq_to_gridxml(imt).lower()
        layer = layers[lname]
        gdata = layer.getData()
        assert np.allclose(gdata, cdata)

    #
    # Do the uncertainty grids
    #
    uxml = os.path.join(datapath, evid, 'current', 'products',
                        'uncertainty.xml')
    u2d = ShakeGrid.load(uxml)
    ulayers = u2d.getData()

    for imt in imts:
        component, imt = imt.split('/')
        comp = oc.getComponents(imt)
        cdata = oc.getIMTGrids(imt, comp[0])['std']
        #
        # The stddevs just get rounded
        #
        digits = oc.getIMTGrids(imt, comp[0])['std_metadata']['digits']
        vfunc = rounder(digits)
        cdata = vfunc(cdata)
        lname = 'std' + _oq_to_gridxml(imt).lower()
        layer = ulayers[lname]
        gdata = layer.getData()
        assert np.allclose(gdata, cdata)


def do_contour(evid, datapath):

    check_failures(evid, datapath, ContourModule)

    mod = ContourModule(evid)
    mod.execute()
    mod.writeContents()

    tfiles = glob.glob(os.path.join(datapath, '..', '..', 'shakemap',
                                    'coremods', 'data', evid, 'cont_*.json'))
    for ifile in tfiles:
        with open(ifile, 'r') as if1:
            cjson = json.load(if1)
        fname = os.path.basename(ifile)
        jfile = os.path.join(datapath, evid, 'current', 'products', fname)
        with open(jfile, 'r') as if2:
            fjson = json.load(if2)
        assert ordered(fjson) == ordered(cjson)


def do_contour_command_line(evid, datapath):
    #
    # Now run from the command line to exercise the argument parsing
    #
    cp = subprocess.run(['shake', '--force', evid, 'contour', '-f', '10'],
                        shell=False)
    assert not cp.returncode


def test_points():
    installpath, datapath = get_config_paths()

    # An event with points rather than grid (should raise exception)
    evid = 'northridge_points'
    assemble = AssembleModule(evid, comment='Test comment.')
    assemble.execute()
    del assemble
    model = ModelModule(evid)
    model.execute()
    del model

    mod = MappingModule(evid)
    # Mapping should raise exception
    with pytest.raises(NotImplementedError):
        mod.execute()


def test_products():

    installpath, datapath = get_config_paths()

    #
    # Use a real event for checks of products against the contents of
    # the output container
    #
    evid = 'integration_test_0001'
    try:
        #
        # Make sure an output file exists
        #
        assemble = AssembleModule(evid, comment='Test comment.')
        assemble.execute()
        del assemble
        model = ModelModule(evid)
        model.execute()
        del model

        res_file = os.path.join(datapath, evid, 'current', 'products',
                                'shake_result.hdf')
        oc = ShakeMapOutputContainer.load(res_file)

        #
        # The history module just outputs some info to the operator, so
        # here we just run it to make sure it doesn't crash anything.
        # Actual testing should be done via bug reports/feature requests
        # from users.
        #
        history = HistoryModule(evid)
        history.execute()
        del history

        #
        # Test the creation of products -- currently not checking results
        # for validity or consistency, but probably should
        #

        #
        # TODO: The stationlist.json should be validated, but we need a
        # function that will read it and convert it to something
        # we can test against.
        #
        check_failures(evid, datapath, StationModule)
        mod = StationModule(evid)
        mod.execute()
        mod.writeContents()

        check_failures(evid, datapath, MappingModule)
        mod = MappingModule(evid)
        mod.execute()
        mod.writeContents()

        check_failures(evid, datapath, PlotRegr)
        #
        # PlotRegr gets tested in the model tests for event 72282711
        #
#        mod = PlotRegr(evid)
#        mod.execute()
#        mod.writeContents()

        check_failures(evid, datapath, KMLModule)
        mod = KMLModule(evid)
        mod.execute()
        mod.writeContents()
        del mod

        # This just exercises the ShapeModule code without actually
        # checking for valid results.
        mod = ShapeModule(evid)
        mod.execute()
        mod.writeContents()
        del mod

        #
        # These check that the results are consistent with the output
        # container
        #
        do_rupture(evid, datapath, oc)

        # do_info(evid, datapath, oc)

        do_raster(evid, datapath, oc)

        do_gridxml(evid, datapath, oc)

        oc.close()
        #
        # Checks contours against saved versions; if something
        # changes, will need to update the files in
        # data/integration_test_0001
        #
        do_contour(evid, datapath)
#        do_contour_command_line(evid, datapath)

    finally:
        pass
        data_file = os.path.join(datapath, evid, 'current', 'shake_data.hdf')
        if os.path.isfile(data_file):
            os.remove(data_file)


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_products()
    test_points()
