#!/usr/bin/env python

# stdlib modules
import sys
import os.path
import time as time
import shutil

# third party modules
import numpy as np
from openquake.hazardlib.gsim import base
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
import openquake.hazardlib.const as oqconst
from openquake.hazardlib.imt import MMI, PGA, PGV, SA
import pytest

# hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
# put this at the front of the system path, ignoring any installed
# shakemap stuff
sys.path.insert(0, shakedir)

# local imports
from shakemap.grind.virtualipe import VirtualIPE
from shakemap.grind.gmice.wgrw12 import WGRW12
from shakemap.grind.multigmpe import MultiGMPE
from shakemap.grind.distance import Distance
from shakemap.grind.sites import Sites
from shakemap.grind.source import Source
from shakemap.utils.exception import ShakeMapException

def test_virtualipe():

    #
    # Set up the GMPE, IPE, and GMICE
    #
    gmpe_cy14 = ChiouYoungs2014()
    gmpe = MultiGMPE.from_list([gmpe_cy14], [1.0])
    gmice = WGRW12()
    ipe = VirtualIPE.fromFuncs(gmpe, gmice)

    #
    # Use the Calexico event info
    #
    homedir = os.path.dirname(os.path.abspath(__file__))
    datadir = os.path.abspath(os.path.join(homedir, '..', 'data',
            'eventdata', 'Calexico', 'input'))

    #
    # Read the event, source, and fault files and produce a Source object
    #
    inputfile = os.path.join(datadir, 'stationlist_dat.xml')
    dyfifile = os.path.join(datadir, 'ciim3_dat.xml')
    eventfile = os.path.join(datadir, 'event.xml')
    faultfile = os.path.join(datadir, 'wei_fault.txt')

    source_obj = Source.fromFile(eventfile, faultfile=faultfile)
    rx = source_obj.getRuptureContext([gmpe])

    smdx = 0.0083333333
    smdy = 0.0083333333
    lonspan = 6.0
    latspan = 4.0
    vs30filename = os.path.join(datadir, '..', 'vs30', 'vs30.grd')

    sites_obj_grid = Sites.fromCenter(
            rx.hypo_lon, rx.hypo_lat, lonspan, latspan,
            smdx, smdy, defaultVs30=760.0, vs30File=vs30filename,
            vs30measured_grid=None, padding=False, resample=False
        )

    npts = 200
    lats = np.empty(npts)
    lons = np.empty(npts)
    depths = np.zeros(npts)
    for i in range(npts):
        lats[i] = rx.hypo_lat
        lons[i] = rx.hypo_lon + i * 0.01
    lldict = {'lats': lats, 'lons': lons}

    sx = sites_obj_grid.getSitesContext(lldict=lldict, rock_vs30=760.0)

    dobj = Distance(gmpe, source_obj, lats, lons, depths)
    dx = dobj.getDistanceContext()

    sd_types = [oqconst.StdDev.TOTAL]
    mmi_const_vs30, mmi_sd_const_vs30 = \
            ipe.get_mean_and_stddevs(sx, rx, dx, MMI(), sd_types)

# These prints are just so a human can examine the outputs
#    print(mmi_const_vs30)
#    print(mmi_sd_const_vs30)

    sx = sites_obj_grid.getSitesContext(lldict=lldict)
    mmi_variable_vs30, mmi_sd_variable_vs30 = \
            ipe.get_mean_and_stddevs(sx, rx, dx, MMI(), sd_types)

#    print(mmi_variable_vs30)
#    print(mmi_sd_variable_vs30)

    sd_types = [oqconst.StdDev.INTRA_EVENT]
    mmi_variable_vs30_intra, mmi_sd_variable_vs30_intra = \
            ipe.get_mean_and_stddevs(sx, rx, dx, MMI(), sd_types)

#    print(mmi_variable_vs30_intra)
#    print(mmi_sd_variable_vs30_intra)
#    assert(0)      # Assert causes test to fail and prints to be displayed

    #
    # Try with PGA
    #
    gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES.remove(PGV)
    ipe = VirtualIPE.fromFuncs(gmpe, gmice)
    mmi_pga, mmi_sd_pga = \
            ipe.get_mean_and_stddevs(sx, rx, dx, MMI(), sd_types)
    #
    # Try with SA(1.0)
    #
    gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES.remove(PGA)
    ipe = VirtualIPE.fromFuncs(gmpe, gmice)
    mmi_psa, mmi_sd_psa = \
            ipe.get_mean_and_stddevs(sx, rx, dx, MMI(), sd_types)

    #
    # This should raise an exception because the IMT isn't MMI
    #
    with pytest.raises(ValueError) as e:
        mmi_psa, mmi_sd_psa = \
                ipe.get_mean_and_stddevs(sx, rx, dx, PGA(), sd_types)
    #
    # This should raise an exception because no valid IMTs are available
    #
    gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES.remove(SA)
    with pytest.raises(ShakeMapException) as e:
        ipe = VirtualIPE.fromFuncs(gmpe, gmice)

    #
    # Now do a GMPE that uses Rjb instead of Rrup
    #
    gmpe_ba14 = BooreEtAl2014()
    gmpe = MultiGMPE.from_list([gmpe_ba14], [1.0])
    ipe = VirtualIPE.fromFuncs(gmpe, gmice)
    rx = source_obj.getRuptureContext([gmpe])
    dobj = Distance(gmpe, source_obj, lats, lons, depths)
    dx = dobj.getDistanceContext()

    mmi_rjb, mmi_sd_rjb = \
            ipe.get_mean_and_stddevs(sx, rx, dx, MMI(), sd_types)
    
    #
    # Test the results against a known standard
    #
    savefile = os.path.abspath(os.path.join(homedir, '..', 'data',
            'eventdata', 'Calexico', 'virtualipe_test', 'savefile.npz'))

    #
    # If things change, set remake_save to True, and it will rebuild the
    # saved data file against which the comparisons are done
    # Remember to set this back to False once you've remade the test datafile
    #
    remake_save = False
    if remake_save:
        np.savez_compressed(savefile,
                mmi_const_vs30 = mmi_const_vs30,
                mmi_sd_const_vs30 = mmi_sd_const_vs30,
                mmi_variable_vs30 = mmi_variable_vs30,
                mmi_sd_variable_vs30 = mmi_sd_variable_vs30,
                mmi_variable_vs30_intra = mmi_variable_vs30_intra,
                mmi_sd_variable_vs30_intra = mmi_sd_variable_vs30_intra,
                mmi_pga = mmi_pga,
                mmi_sd_pga = mmi_sd_pga,
                mmi_psa = mmi_psa,
                mmi_sd_psa = mmi_sd_psa,
                mmi_rjb = mmi_rjb,
                mmi_sd_rjb = mmi_sd_rjb)

    td = np.load(savefile)

    assert(np.allclose(td['mmi_const_vs30'], mmi_const_vs30))
    assert(np.allclose(td['mmi_sd_const_vs30'], mmi_sd_const_vs30))
    assert(np.allclose(td['mmi_variable_vs30'], mmi_variable_vs30))
    assert(np.allclose(td['mmi_sd_variable_vs30'], mmi_sd_variable_vs30))
    assert(np.allclose(td['mmi_variable_vs30_intra'], mmi_variable_vs30_intra))
    assert(np.allclose(td['mmi_sd_variable_vs30_intra'], 
        mmi_sd_variable_vs30_intra))
    assert(np.allclose(td['mmi_pga'], mmi_pga))
    assert(np.allclose(td['mmi_sd_pga'], mmi_sd_pga))
    assert(np.allclose(td['mmi_psa'], mmi_psa))
    assert(np.allclose(td['mmi_sd_psa'], mmi_sd_psa))
    assert(np.allclose(td['mmi_rjb'], mmi_rjb))
    assert(np.allclose(td['mmi_sd_rjb'], mmi_sd_rjb))

    # The total uncertainties should be greater than the intra-event
    assert(np.all(mmi_sd_variable_vs30 > mmi_sd_variable_vs30_intra))
    
if __name__ == '__main__':

    test_virtualipe()
