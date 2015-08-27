#!/usr/bin/env python

from openquake.hazardlib.gsim.boore_atkinson_2008 import BooreAtkinson2008
from openquake.hazardlib import const

import numpy as np

class BooreAtkinson2008ShakeMap(BooreAtkinson2008):
    '''
    Subclasses GEM OpenQuake Boore-Atkinson 2008 GMPE, but with methods that allow finer-grained
    access to methods into that GMPE implementation.  Specifically, methods that allow the user to:
     - apply or un-apply site corrections to a particular site.
     - get ground motions on bare rock
     - get the standard deviations for the GMPE
    '''

    def get_site_corrections(self, sites, rup, dists, imt, pgm, forward=True):
        '''
        Calculate site corrections for sites on rock (forward=True), OR
        Remove site corrections from data on other substrates.
        '''
        C = self.COEFFS[imt]
        C_SR = self.COEFFS_SOIL_RESPONSE[imt]
        pga4nl = self._get_pga_on_rock(rup, dists, C)
        if forward:
            pgm_corrected = pgm + \
               self._get_site_amplification_linear(sites.vs30, C_SR) + \
               self._get_site_amplification_non_linear(sites.vs30, pga4nl, C_SR)
        else:
            pgm_corrected = pgm - \
               (self._get_site_amplification_linear(sites.vs30, C_SR) + \
                self._get_site_amplification_non_linear(sites.vs30, pga4nl, C_SR))

        return pgm_corrected

    def get_amplitudes(self, rup, dists, imt):
        '''
        Calculate peak ground motion on rock.
        '''
        C = self.COEFFS[imt]
        pgmrock = self._compute_magnitude_scaling(rup, C) + \
                  self._compute_distance_scaling(rup, dists, C)

        return pgmrock

    def get_stddevs(self, sites, rup, dists, imt, stddev_type):
        '''
        Get standard deviations of GMPE.
        '''
        C = self.COEFFS[imt]

        assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
        if stddev_type == const.StdDev.TOTAL:
            stddevs = np.zeros_like(dists.rjb) + C['std']
        elif stddev_type == const.StdDev.INTRA_EVENT:
            stddevs = np.zeros_like(dists.rjb) + C['sigma']
        elif stddev_type == const.StdDev.INTER_EVENT:
            stddevs = np.zeros_like(dists.rjb) + C['tau']
        else:
            raise ValueError("Unknown type %s" % stddev_type)
        return stddevs

    def _get_fault_type_dummy_variables(self, rup):
        """
        Override this function to allow for undefined rupture type
        """
        if rup.rake is None:
            return 1, 0, 0, 0
        return super(BooreAtkinson2008ShakeMap, self)._get_fault_type_dummy_variables(rup)
