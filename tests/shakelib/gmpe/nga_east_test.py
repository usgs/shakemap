#!/usr/bin/env python

import os
import pickle

import numpy as np

from openquake.hazardlib.gsim import base
import openquake.hazardlib.imt as imt
from openquake.hazardlib.const import StdDev

from shakelib.gmpe.nga_east import NGAEast

home_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(home_dir, 'nga_east_data')

stddev_types = [StdDev.TOTAL]
gmpe = NGAEast()

dx = base.DistancesContext()
dx.rrup = np.logspace(-1, np.log10(2000), 100)

rx = base.RuptureContext()
sx = base.SitesContext()

IMTS = [
    imt.PGA(),
    imt.PGV(),
    imt.SA(0.3),
    imt.SA(1.0),
    imt.SA(3.0)
]

MAGS = [3, 5, 6, 7]

VS30 = [180, 380, 760, 2000]


def update_results():
    # To build the data for testing
    result = {}
    for i in IMTS:
        ikey = i.__str__()
        result[ikey] = {}
        for mag in MAGS:
            rx.mag = mag
            result[ikey][str(mag)] = {}
            for vs30 in VS30:
                sx.vs30 = np.full_like(dx.rrup, vs30)
                result[ikey][str(mag)][str(vs30)] = {}
                lmean, lsd = gmpe.get_mean_and_stddevs(
                    sx, rx, dx, i, stddev_types)
                result[ikey][str(mag)][str(vs30)]['lmean'] = lmean.tolist()
                result[ikey][str(mag)][str(vs30)]['lsd'] = lsd[0].tolist()
    # Save results
    pkl_file = os.path.join(data_dir, 'nga_east_data.pkl')
    fh = open(pkl_file, 'wb')
    pickle.dump(result, fh)
    fh.close()


def test_nga_east():
    # Load test data
    pkl_file = os.path.join(data_dir, 'nga_east_data.pkl')
    fh = open(pkl_file, 'rb')
    target = pickle.load(fh)
    fh.close()
    for i in IMTS:
        ikey = i.__str__()
        for mag in MAGS:
            rx.mag = mag
            for vs30 in VS30:
                sx.vs30 = np.full_like(dx.rrup, vs30)
                lmean, lsd = gmpe.get_mean_and_stddevs(
                    sx, rx, dx, i, stddev_types)
                tmean = np.array(target[ikey][str(mag)][str(vs30)]['lmean'])
                np.testing.assert_allclose(lmean, tmean, rtol=1e-6, atol=1e-6)
                tsd = np.array(target[ikey][str(mag)][str(vs30)]['lsd'])
                np.testing.assert_allclose(lsd[0], tsd, rtol=1e-6, atol=1e-6)


if __name__ == '__main__':
    test_nga_east()
