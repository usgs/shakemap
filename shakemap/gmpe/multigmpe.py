#!/usr/bin/env python

import numpy as np

from openquake.hazardlib.gsim.base import GMPE
from openquake.hazardlib import const

from shakemap.grind.BeyerBommer2006 import ampIMCtoIMC, sigmaIMCtoIMC


class MultiGMPE(GMPE):
    """Implements a GMPE that is the combination of multiple GMPEs.
    TODO:

       * convert IMT (e.g., PGV) from another IMT if it is not available from
          the GMPE in get_mean_and_stddevs.

    """

    DEFINED_FOR_TECTONIC_REGION_TYPE = None
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = None
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = None
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = None
    REQUIRES_SITES_PARAMETERS = None
    REQUIRES_RUPTURE_PARAMETERS = None
    REQUIRES_DISTANCES = None

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        lnmu = np.zeros_like(sites.vs30)
        lnsd2 = np.zeros_like(sites.vs30)
        lmean = [None] * len(self.GMPEs)
        lsd = [None] * len(self.GMPEs)
        for i in range(len(self.GMPEs)):
            gmpe = self.GMPEs[i]
            # Need to select the appropriate z1pt0 value for different GMPEs.
            # Note that these are required site parameters, so even though
            # OQ has these equations built into the class, the arrays must
            # be provided in the sites context. It might be worth sending
            # a request to OQ to provide a subclass that that computes the
            # depth parameters when not provided (as is done for BSSA14 but
            # not the others).
            if gmpe == 'AbrahamsonEtAl2014()':
                sites.z1pt0 = sites.z1pt0ask14
            if gmpe == 'BooreEtAl2014()' or gmpe == 'ChiouYoungs2014()':
                sites.z1pt0 = sites.z1pt0cy14
            lmean[i], lsd[i] = gmpe.get_mean_and_stddevs(
                sites, rup, dists, imt, stddev_types)

            # Convert component type.
            # Note: conversion is based on linear amps (not log)!!
            inc_in = self.IMCs[i]
            inc_out = self.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT
            lmean[i] = np.log(
                ampIMCtoIMC(
                    np.exp(
                        lmean[i]),
                    inc_in,
                    inc_out,
                    imt))
            lsd[i] = np.log(
                sigmaIMCtoIMC(
                    np.exp(
                        lsd[i]),
                    inc_in,
                    inc_out,
                    imt))

            # Compute weighted mean and sd
            lnmu = lnmu + self.weights[i] * lmean[i]
            lnsd2 = lnsd2 + self.weights[i] * (lmean[i]**2 + lsd[i]**2)
        lnsd2 = lnsd2 - lnmu**2

        return lnmu, np.sqrt(lnsd2)

    @classmethod
    def from_list(cls, GMPEs, weights):
        """Construct a MultiGMPE from lists of GMPEs and weights.

        :param GMPEs:
            List of OpenQuake GMPE instances.
        :param weights:
            List of weights.
        """
        # Check that weights sum to 1.0:
        if np.sum(weights) != 1.0:
            raise Exception('Weights must sum to one.')

        # Check that GMPEs is a list of OQ GMPE instances
        for g in GMPEs:
            if not isinstance(g, GMPE):
                raise Exception("\"%s\" is not a GMPE instance." % g)

        self = cls()
        self.GMPEs = GMPEs
        self.weights = weights

        # Check that GMPEs all are for the same tectonic region,
        # otherwise raise exception.
        tmp = set([i.DEFINED_FOR_TECTONIC_REGION_TYPE for i in GMPEs])
        if len(tmp) == 1:
            self.DEFINED_FOR_TECTONIC_REGION_TYPE = \
                GMPEs[0].DEFINED_FOR_TECTONIC_REGION_TYPE
        else:
            raise Exception('GMPEs are not all for the same tectonic region.')

        # Combine the intensity measure types. This is problematic:
        #   - Logically, we should only include the intersection of the sets
        #     of imts for the different GMPEs.
        #   - In practice, this is not feasible because most GMPEs in CEUS and
        #     subduction zones do not have PGV.
        #   - So instead we will use the union of the imts and then convert
        #     to get the missing imts later in get_mean_and_stddevs.
        imts = [g.DEFINED_FOR_INTENSITY_MEASURE_TYPES for g in GMPEs]
        self.DEFINED_FOR_INTENSITY_MEASURE_TYPES = set.union(*imts)

        # Store intensity measure types for conversion in get_mean_and_stddevs.
        self.IMCs = [g.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT for g in GMPEs]

        # For ShakeMap, the target IMC is max
        self.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = \
            const.IMC.GREATER_OF_TWO_HORIZONTAL

        # For scenarios, we only care about total standard deviation,
        # but for real-time we need inter and intra. For now, lets
        # just take the intersection of the different GMPEs to make life
        # slightly easier.
        stdlist = [set(g.DEFINED_FOR_STANDARD_DEVIATION_TYPES) for g in GMPEs]
        self.DEFINED_FOR_STANDARD_DEVIATION_TYPES = set.intersection(*stdlist)

        # Need union of site parameters, but it is complicated by the
        # different depth parameter flavors.
        sitepars = [g.REQUIRES_SITES_PARAMETERS for g in GMPEs]
        self.REQUIRES_SITES_PARAMETERS = set.union(*sitepars)

        # Union of rupture parameters
        ruppars = [g.REQUIRES_RUPTURE_PARAMETERS for g in GMPEs]
        self.REQUIRES_RUPTURE_PARAMETERS = set.union(*ruppars)

        # Union of distance parameters
        distpars = [g.REQUIRES_DISTANCES for g in GMPEs]
        self.REQUIRES_DISTANCES = set.union(*distpars)

        return self

#----------------------------------------------------
# Functions for getting depth parameters from Vs30
#----------------------------------------------------


def z1_from_vs30_cy14_cal(vs30):
    # vs30 = V_S30 in units of m/s
    # z1   = z_1 in units of m
    z1 = np.exp(-(7.15 / 4.0) *
                np.log((vs30**4.0 + 571.**4) / (1360**4.0 + 571.**4)))
    return z1


def z1_from_vs30_ask14_cal(vs30):
    # vs30 = V_S30 in units of m/s
    # z1   = z_1 in units of m
    # ASK14 define units as KM, but implemented as m in OQ
    z1 = np.exp(-(7.67 / 4.0) *
                np.log((vs30**4.0 + 610.**4) / (1360**4.0 + 610.**4)))
    return z1


def z2p5_from_vs30_cb14_cal(vs30):
    # vs30 = V_S30 in units of m/s
    # z2p5 = z_2.5 in units of m
    z2p5 = 1000 * np.exp(7.089 - (1.144) * np.log(vs30))
    return z2p5


