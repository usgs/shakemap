#!/usr/bin/env python

# stdlib modules
import sys
import os.path

# third party modules
import numpy as np
from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
import openquake.hazardlib.const as oqconst
from openquake.hazardlib.imt import MMI, PGA, PGV, SA
import pytest

# local imports
from shakelib.distance import Distance
from shakelib.gmice.wgrw12 import WGRW12
from shakelib.multigmpe import MultiGMPE
from shakelib.rupture.origin import Origin
from shakelib.rupture.factory import get_rupture
from shakelib.sites import Sites
from shakelib.utils.exception import ShakeLibException
from shakelib.virtualipe import VirtualIPE


homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)


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
    datadir = os.path.abspath(os.path.join(
        homedir, 'virtualipe_data', 'Calexico', 'input'))

    #
    # Read the event, origin, and rupture files and produce Rupture and Origin
    # objects
    #
    # inputfile = os.path.join(datadir, 'stationlist_dat.xml')
    # dyfifile = os.path.join(datadir, 'ciim3_dat.xml')
    eventfile = os.path.join(datadir, 'event.xml')
    rupturefile = os.path.join(datadir, 'wei_fault.txt')

    origin_obj = Origin.fromFile(eventfile)
    rupture_obj = get_rupture(origin_obj, rupturefile)
    rx = rupture_obj.getRuptureContext([gmpe])
    rx.rake = 45.

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

    dobj = Distance(gmpe, lons, lats, depths, rupture_obj)
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

    sd_types = [oqconst.StdDev.TOTAL, oqconst.StdDev.INTRA_EVENT,
                oqconst.StdDev.INTER_EVENT]
    mmi_variable_vs30_intra, mmi_sd_variable_vs30_intra = \
        ipe.get_mean_and_stddevs(sx, rx, dx, MMI(), sd_types)

#    print(mmi_variable_vs30_intra)
#    print(mmi_sd_variable_vs30_intra)
#    assert(0)      # Assert causes test to fail and prints to be displayed

    #
    # Try with PGA
    #
    gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES.remove(PGV)
    gmpe.ALL_GMPES_HAVE_PGV = False
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
    with pytest.raises(ShakeLibException) as e:  # noqa
        ipe = VirtualIPE.fromFuncs(gmpe, gmice)

    #
    # Now do a GMPE that uses Rjb instead of Rrup
    #
    gmpe_ba14 = BooreEtAl2014()
    gmpe = MultiGMPE.from_list([gmpe_ba14], [1.0])
    ipe = VirtualIPE.fromFuncs(gmpe, gmice)
    rx = rupture_obj.getRuptureContext([gmpe])
    rx.rake = 45.
    dobj = Distance(gmpe, lons, lats, depths, rupture_obj)
    dx = dobj.getDistanceContext()

    mmi_rjb, mmi_sd_rjb = \
        ipe.get_mean_and_stddevs(sx, rx, dx, MMI(), sd_types)

    #
    # Test the results against a known standard
    #
    savefile = os.path.abspath(os.path.join(
        homedir, 'virtualipe_data', 'Calexico', 'virtualipe_test',
        'savefile.npz'))

    #
    # If things change, set remake_save to True, and it will rebuild the
    # saved data file against which the comparisons are done
    # Remember to set this back to False once you've remade the test datafile
    #
    remake_save = False
    if remake_save:
        np.savez_compressed(
            savefile,
            mmi_const_vs30=mmi_const_vs30,
            mmi_sd_const_vs30=mmi_sd_const_vs30[0],
            mmi_variable_vs30=mmi_variable_vs30,
            mmi_sd_variable_vs30=mmi_sd_variable_vs30[0],
            mmi_variable_vs30_intra=mmi_variable_vs30_intra,
            mmi_sd_variable_vs30_total=mmi_sd_variable_vs30_intra[0],
            mmi_sd_variable_vs30_intra=mmi_sd_variable_vs30_intra[1],
            mmi_sd_variable_vs30_inter=mmi_sd_variable_vs30_intra[2],
            mmi_pga=mmi_pga,
            mmi_sd_pga=mmi_sd_pga[0],
            mmi_psa=mmi_psa,
            mmi_sd_psa=mmi_sd_psa[0],
            mmi_rjb=mmi_rjb,
            mmi_sd_rjb=mmi_sd_rjb[0])

    td = np.load(savefile)

    assert(np.allclose(td['mmi_const_vs30'], mmi_const_vs30))
    assert(np.allclose(td['mmi_sd_const_vs30'], mmi_sd_const_vs30[0]))
    assert(np.allclose(td['mmi_variable_vs30'], mmi_variable_vs30))
    assert(np.allclose(td['mmi_sd_variable_vs30'], mmi_sd_variable_vs30[0]))
    assert(np.allclose(td['mmi_variable_vs30_intra'], mmi_variable_vs30_intra))
    assert(np.allclose(td['mmi_sd_variable_vs30_total'],
                       mmi_sd_variable_vs30_intra[0]))
    assert(np.allclose(td['mmi_sd_variable_vs30_intra'],
                       mmi_sd_variable_vs30_intra[1]))
    assert(np.allclose(td['mmi_sd_variable_vs30_inter'],
                       mmi_sd_variable_vs30_intra[2]))
    assert(np.allclose(td['mmi_pga'], mmi_pga))
    assert(np.allclose(td['mmi_sd_pga'], mmi_sd_pga[0]))
    assert(np.allclose(td['mmi_psa'], mmi_psa))
    assert(np.allclose(td['mmi_sd_psa'], mmi_sd_psa[0]))
    assert(np.allclose(td['mmi_rjb'], mmi_rjb))
    assert(np.allclose(td['mmi_sd_rjb'], mmi_sd_rjb[0]))

    # The total uncertainties should be greater than the intra-event
    assert(np.all(mmi_sd_variable_vs30[0] > mmi_sd_variable_vs30_intra[1]))

    # The combined intra and inter-event uncertainty should be equal
    # to the total
    tot = np.sqrt(
        mmi_sd_variable_vs30_intra[1]**2 + mmi_sd_variable_vs30_intra[2]**2)
    assert(np.allclose(tot, mmi_sd_variable_vs30_intra[0], rtol=1e-2))

    # Run through mgm = mgm + fd if fd is available
    # In this case fd = 0, so the ouput is still the same values
    mmi_fd, mmi_fdsd = ipe.get_mean_and_stddevs(sx, rx, dx, MMI(),
                                                sd_types, fd=0)
    assert(np.allclose(mmi_fd, mmi_rjb))
    assert(np.allclose(mmi_fdsd, mmi_sd_rjb))


if __name__ == '__main__':

    test_virtualipe()
