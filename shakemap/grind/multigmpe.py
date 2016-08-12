#!/usr/bin/env python

import numpy as np

from openquake.hazardlib.gsim.base import GMPE
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGV
from openquake.hazardlib.imt import SA

from shakemap.grind.conversions.imt.newmark_hall_1982 import NewmarkHall1982
from shakemap.grind.conversions.imc.beyer_bommer_2006 import BeyerBommer2006


class MultiGMPE(GMPE):
    """
    Implements a GMPE that is the combination of multiple GMPEs.

    To do

        * Update IMT conversion to account for additional uncertainty. 
        * Develop a method to include GMPEs that don't have a site term. 
        * Add a check that the lenght of the weights match the length 
          of the GMPE list. 

    """

    DEFINED_FOR_TECTONIC_REGION_TYPE = None
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = None
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = None
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = None
    REQUIRES_SITES_PARAMETERS = None
    REQUIRES_RUPTURE_PARAMETERS = None
    REQUIRES_DISTANCES = None

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See superclass `method <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#openquake.hazardlib.gsim.base.GroundShakingIntensityModel.get_mean_and_stddevs>`__. 
        """

        # These are arrays to hold the weighted combination of the GMPEs
        lnmu = np.zeros_like(sites.vs30)
        stddev_types = self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
        lnsd2 = [np.zeros_like(sites.vs30) for a in stddev_types]

        for i in range(len(self.GMPEs)):
            #---------------------------------------------------------------
            # Loop over GMPE list
            #---------------------------------------------------------------

            gmpe = self.GMPEs[i]

            #---------------------------------------------------------------
            # Need to select the appropriate z1pt0 value for different GMPEs.
            # Note that these are required site parameters, so even though
            # OQ has these equations built into the class, the arrays must
            # be provided in the sites context. It might be worth sending
            # a request to OQ to provide a subclass that that computes the
            # depth parameters when not provided (as is done for BSSA14 but
            # not the others).
            #---------------------------------------------------------------

            if gmpe == 'AbrahamsonEtAl2014()':
                sites.z1pt0 = sites.z1pt0ask14
            if gmpe == 'BooreEtAl2014()' or gmpe == 'ChiouYoungs2014()':
                sites.z1pt0 = sites.z1pt0cy14


            #---------------------------------------------------------------
            # Evaluate GMPEs
            #---------------------------------------------------------------

            gmpe_imts = [imt.__name__ for imt in \
                         gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES]
            if (isinstance(imt, PGV)) and ("PGV" not in gmpe_imts):

                # If IMT is PGV and not given by GMPE, convert from PSA10

                psa10, psa10sd = gmpe.get_mean_and_stddevs(
                    sites, rup, dists, SA(1.0), stddev_types)
                lmean, lsd = NewmarkHall1982.psa102pgv(psa10, psa10sd[0])
            else:
                lmean, lsd = gmpe.get_mean_and_stddevs(
                    sites, rup, dists, imt, stddev_types)

            #---------------------------------------------------------------
            # Convertions due to component definition
            #---------------------------------------------------------------

            imc_in = gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT
            imc_out = self.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT
            lmean = BeyerBommer2006.ampIMCtoIMC(lmean, imc_in, imc_out, imt)
            for j in range(len(lnsd2)):
                lsd[j] = BeyerBommer2006.sigmaIMCtoIMC(
                    lsd[j], imc_in, imc_out, imt)

            #---------------------------------------------------------------
            # Compute weighted mean and sd
            #---------------------------------------------------------------

            lnmu = lnmu + self.weights[i] * lmean

            # Note: the lnsd2 calculation isn't complete until we drop out of
            # this loop and substract lnmu**2
            for j in range(len(lnsd2)):
                lnsd2[j] = lnsd2[j] + self.weights[i] * (lmean**2 + lsd[j]**2)

        for j in range(len(lnsd2)):
            lnsd2[j] = lnsd2[j] - lnmu**2

        lnsd = [np.sqrt(a) for a in lnsd2]

        return lnmu, lnsd

    @classmethod
    def from_list(cls, GMPEs, weights, imc = const.IMC.GREATER_OF_TWO_HORIZONTAL):
        """Construct a MultiGMPE instance from lists of GMPEs and weights.

        :param GMPEs:
            List of OpenQuake 
            `GMPE <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__ 
            instances.
        :param weights:
            List of weights; must sum to 1.0.
        :param imc: Requested intensity measure component. Must be one listed
            `here <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__.
            The amplitudes returned by the GMPEs will be converted to this IMT. 
            Default is 'GREATER_OF_TWO_HORIZONTAL', which is used by ShakeMap. 
            See discussion in `this section <http://usgs.github.io/shakemap/tg_choice_of_parameters.html#use-of-peak-values-rather-than-mean>`__
            of the ShakeMap manual. 
        """

        #---------------------------------------------------------
        # Check that weights sum to 1.0:
        #---------------------------------------------------------

        if np.sum(weights) != 1.0:
            raise Exception('Weights must sum to one.')

        #---------------------------------------------------------
        # Check that length of weights equals length of gmpe list
        #---------------------------------------------------------
        if len(weights) != len(GMPEs):
            raise Exception('Length of weights must match length of GMPE list.')

        #---------------------------------------------------------
        # Check that GMPEs is a list of OQ GMPE instances
        #---------------------------------------------------------

        for g in GMPEs:
            if not isinstance(g, GMPE):
                raise Exception("\"%s\" is not a GMPE instance." % g)

        self = cls()
        self.GMPEs = GMPEs
        self.weights = weights

        #---------------------------------------------------------
        # Check that GMPEs all are for the same tectonic region,
        # otherwise raise exception.
        #---------------------------------------------------------

        tmp = set([i.DEFINED_FOR_TECTONIC_REGION_TYPE for i in GMPEs])
        if len(tmp) == 1:
            self.DEFINED_FOR_TECTONIC_REGION_TYPE = \
                GMPEs[0].DEFINED_FOR_TECTONIC_REGION_TYPE
        else:
            raise Exception('GMPEs are not all for the same tectonic region.')

        #---------------------------------------------------------
        # Combine the intensity measure types. This is problematic:
        #   - Logically, we should only include the intersection of the sets
        #     of imts for the different GMPEs.
        #   - In practice, this is not feasible because most GMPEs in CEUS and
        #     subduction zones do not have PGV.
        #   - So instead we will use the union of the imts and then convert
        #     to get the missing imts later in get_mean_and_stddevs.
        #---------------------------------------------------------

        imts = [g.DEFINED_FOR_INTENSITY_MEASURE_TYPES for g in GMPEs]
        self.DEFINED_FOR_INTENSITY_MEASURE_TYPES = set.union(*imts)

        #---------------------------------------------------------
        # Store intensity measure types for conversion in get_mean_and_stddevs.
        #---------------------------------------------------------
        self.IMCs = [g.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT for g in GMPEs]

        #---------------------------------------------------------
        # Store the component
        #---------------------------------------------------------
        self.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = imc

        #---------------------------------------------------------
        # For scenarios, we only care about total standard deviation,
        # but for real-time we need inter and intra. For now, lets
        # just use total to make life easier. 
        #---------------------------------------------------------
#        stdlist = [set(g.DEFINED_FOR_STANDARD_DEVIATION_TYPES) for g in GMPEs]
        self.DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
            const.StdDev.TOTAL
        ])

        #---------------------------------------------------------
        # Need union of site parameters, but it is complicated by the
        # different depth parameter flavors.
        #---------------------------------------------------------
        sitepars = [g.REQUIRES_SITES_PARAMETERS for g in GMPEs]
        self.REQUIRES_SITES_PARAMETERS = set.union(*sitepars)

        #---------------------------------------------------------
        # Union of rupture parameters
        #---------------------------------------------------------
        ruppars = [g.REQUIRES_RUPTURE_PARAMETERS for g in GMPEs]
        self.REQUIRES_RUPTURE_PARAMETERS = set.union(*ruppars)

        #---------------------------------------------------------
        # Union of distance parameters
        #---------------------------------------------------------
        distpars = [g.REQUIRES_DISTANCES for g in GMPEs]
        self.REQUIRES_DISTANCES = set.union(*distpars)

        return self

#----------------------------------------------------
# Functions for getting depth parameters from Vs30
#----------------------------------------------------


def _z1_from_vs30_cy14_cal(vs30):
    """
    Compute z1.0 using CY14 relationship. 

    :param vs30:
        Numpy array of Vs30 values in m/s. 
    :returns: 
        Numpy array of z1.0 in m.  
    """
    z1 = np.exp(-(7.15 / 4.0) *
                np.log((vs30**4.0 + 571.**4) / (1360**4.0 + 571.**4)))
    return z1


def _z1_from_vs30_ask14_cal(vs30):
    """
    Calculate z1.0 using ASK14 relationship. 

    :param vs30:
        Numpy array of Vs30 values in m/s. 
    :returns: 
        Numpy array of z1.0 in m.  

    """
    # ASK14 define units as KM, but implemented as m in OQ
    z1 = np.exp(-(7.67 / 4.0) *
                np.log((vs30**4.0 + 610.**4) / (1360**4.0 + 610.**4)))
    return z1


def _z2p5_from_vs30_cb14_cal(vs30):
    """
    Calculate z2.5 using CB14 relationship. 

    :param vs30:
        Numpy array of Vs30 values in m/s. 
    :returns: 
        Numpy array of z2.5 in m.  
    """
    z2p5 = 1000 * np.exp(7.089 - (1.144) * np.log(vs30))
    return z2p5
