"""
Module exports class BB
"""

from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, PGV, SA

import numpy as np

"""
IMC equivalencies:
OpenQuake                   Beyer&Bommer
----------------------------------------------
AVERAGE_HORIZONTAL          GMxy (Geometric mean)       # This is the "reference" type
HORIZONTAL                  ???
MEDIAN_HORIZONTAL           AMxy (Arithmetic mean)
GMRotI50                    GMRotI50
RotD50                      GMRotD50
RANDOM_HORIZONTAL           Random
GREATER_OF_TWO_HORIZONTAL   Env_xy
VERTICAL                    ---
"""

class BB(object):
    """
    Implements conversion for various "Intensity Measure
    Components" per Beyer and Bommer, 2006, BSSA, 96(4A).
    """
    def __init__(self):
        pga_pgv_col_names = [ 'c12', 'c34', 'R' ]
        sa_col_names = [ 'c1', 'c2', 'c3', 'c4', 'R' ]
        self.pga_dict = { \
            const.IMC.MEDIAN_HORIZONTAL         : dict(zip(pga_pgv_col_names, [ 1.0, 0.01, 1.00 ])), \
            const.IMC.GMRotI50                  : dict(zip(pga_pgv_col_names, [ 1.0, 0.02, 1.00 ])), \
            const.IMC.RotD50                    : dict(zip(pga_pgv_col_names, [ 1.0, 0.02, 1.00 ])), \
            const.IMC.RANDOM_HORIZONTAL         : dict(zip(pga_pgv_col_names, [ 1.0, 0.07, 1.03 ])), \
            const.IMC.GREATER_OF_TWO_HORIZONTAL : dict(zip(pga_pgv_col_names, [ 1.1, 0.05, 1.02 ])) }
        self.pgv_dict = { \
            const.IMC.MEDIAN_HORIZONTAL         : dict(zip(pga_pgv_col_names, [ 1.0, 0.01, 1.00 ])), \
            const.IMC.GMRotI50                  : dict(zip(pga_pgv_col_names, [ 1.0, 0.02, 1.00 ])), \
            const.IMC.RotD50                    : dict(zip(pga_pgv_col_names, [ 1.0, 0.02, 1.00 ])), \
            const.IMC.RANDOM_HORIZONTAL         : dict(zip(pga_pgv_col_names, [ 1.0, 0.07, 1.03 ])), \
            const.IMC.GREATER_OF_TWO_HORIZONTAL : dict(zip(pga_pgv_col_names, [ 1.1, 0.05, 1.02 ])) }
        self.sa_dict  = { \
            const.IMC.MEDIAN_HORIZONTAL         : dict(zip(sa_col_names, [ 1.0, 1.0, 0.01, 0.02, 1.00 ])), \
            const.IMC.GMRotI50                  : dict(zip(sa_col_names, [ 1.0, 1.0, 0.03, 0.04, 1.00 ])), \
            const.IMC.RotD50                    : dict(zip(sa_col_names, [ 1.0, 1.0, 0.02, 0.03, 1.00 ])), \
            const.IMC.RANDOM_HORIZONTAL         : dict(zip(sa_col_names, [ 1.0, 1.0, 0.07, 0.11, 1.05 ])), \
            const.IMC.GREATER_OF_TWO_HORIZONTAL : dict(zip(sa_col_names, [ 1.1, 1.2, 0.04, 0.07, 1.02 ])) }

    def ampIMCtoIMC(self, amps, imc_in, imc_out, imt):
        """ 
        Returns amps converted from one IMC to another.
        IMPORTANT: Assumes the input amps are in linear (not log) space 
        IMPORTANT: IMC types 'VERTICAL' and 'HORIZONTAL' are not supported
        Inputs:
            amps -- a numpy array of ground motions in IMC imc_in and IMT imt
            imc_in -- the IMC type of the input amp array
            imc_out -- the desired IMC type of the output amps
            imt -- the IMT of the input amps (must be one of PGA, PGV, or SA)
        Outputs:
            returns amps converted from imc_in to imc_out
        """
        if imc_in == const.IMC.AVERAGE_HORIZONTAL:
            # The amps are already in the B&B "reference" type ("GM", i.e., 
            # geometric mean)
            denom = 1
        elif imc_in == const.IMC.GREATER_OF_TWO_HORIZONTAL or \
             imc_in == const.IMC.MEDIAN_HORIZONTAL or \
             imc_in == const.IMC.GMRotI50          or\
             imc_in == const.IMC.RotD50            or\
             imc_in == const.IMC.RANDOM_HORIZONTAL:
            denom = self._GM2other(imt, imc_in)
        else:
            raise ValueError('unknown IMC %r' % imc_in)

        if imc_out == const.IMC.AVERAGE_HORIZONTAL:
            # The previous step will put the amps into the B&B "reference" 
            # type ("GM", i.e. geometric mean)
            numer = 1
        elif imc_out == const.IMC.GREATER_OF_TWO_HORIZONTAL or \
             imc_out == const.IMC.MEDIAN_HORIZONTAL or \
             imc_out == const.IMC.GMRotI50          or\
             imc_out == const.IMC.RotD50            or\
             imc_out == const.IMC.RANDOM_HORIZONTAL:
            numer = self._GM2other(imt, imc_out)
        else:
            raise ValueError('unknown IMC %r' % imc_out)

        return amps * (numer / denom)

    def sigmaIMCtoIMC(self, sigmas, imc_in, imc_out, imt):
        """ 
        Returns sigmas converted from one IMC to another.
        IMPORTANT: Assumes the input sigmas are in log space 
        IMPORTANT: IMC types 'VERTICAL' and 'HORIZONTAL' are not supported
        Inputs:
            sigmas -- a numpy array of ground motions in IMC imc_in and IMT imt
            imc_in -- the IMC type of the input sigmas array
            imc_out -- the desired IMC type of the output sigmas
            imt -- the IMT of the input sigmas (must be one of PGA, PGV, or SA)
        Outputs:
            returns sigmas converted from imc_in to imc_out
        """
        if imc_in == const.IMC.AVERAGE_HORIZONTAL:
            # The amps are already in the B&B "reference" type ("GM", i.e., 
            # geometric mean)
            R = 1
            sig_log_ratio = 0
        elif imc_in == const.IMC.GREATER_OF_TWO_HORIZONTAL or \
             imc_in == const.IMC.MEDIAN_HORIZONTAL or \
             imc_in == const.IMC.GMRotI50          or\
             imc_in == const.IMC.RotD50            or\
             imc_in == const.IMC.RANDOM_HORIZONTAL:
            R, sig_log_ratio = self._GM2otherSigma(imt, imc_in)
        else:
            raise ValueError('unknown IMC %r' % imc_in)

        sigma_GM2 = (sigmas**2 - sig_log_ratio**2) / R**2

        if imc_out == const.IMC.AVERAGE_HORIZONTAL:
            # The sigmas (sigma_GM2) are already in the B&B "reference" type ("GM", i.e., 
            # geometric mean)
            R = 1
            sig_log_ratio = 0
        elif imc_out == const.IMC.GREATER_OF_TWO_HORIZONTAL or \
             imc_out == const.IMC.MEDIAN_HORIZONTAL or \
             imc_out == const.IMC.GMRotI50          or\
             imc_out == const.IMC.RotD50            or\
             imc_out == const.IMC.RANDOM_HORIZONTAL:
            R, sig_log_ratio = self._GM2otherSigma(imt, imc_out)
        else:
            raise ValueError('unknown IMC %r' % imc_out)

        sigma_out = np.sqrt(sigma_GM2 * R**2 + sig_log_ratio**2)
        return sigma_out

    def _GM2other(self, imt, imc):
        """ Helper function to extract coefficients from the parameter tables """
        if 'PGA' in imt:
            return self.pga_dict[imc]['c12']
        elif 'PGV' in imt:
            return self.pgv_dict[imc]['c12']
        elif 'SA' in imt:
            pp = imt.period
            if pp <= 0.15:
                return self.sa_dict[imc]['c1']
            elif pp < 0.8:
                c1 = self.sa_dict[imc]['c1']
                c2 = self.sa_dict[imc]['c2']
                return c1 + (c2 - c1) * log(pp / 0.15) / log(0.8 / 0.15)
            elif pp <= 5.0:
                return self.sa_dict[imc]['c2']
            else:
                # Not sure what's right here; should probably raise an error
                # but for now let's just use c2
                return self.sa_dict[imc]['c2']
        else:
            raise ValueError('unknown IMT %r' % imt)

    def _GM2otherSigma(self, imt, imc):
        """ Helper function to extract coefficients from the parameter tables """
        if 'PGA' in imt:
            return self.pga_dict[imc]['R'], self.pga_dict[imc]['c34']
        elif 'PGV' in imt:
            return self.pgv_dict[imc]['R'], self.pgv_dict[imc]['c34']
        elif 'SA' in imt:
            R = self.sa_dict[imc]['R']
            pp = imt.period
            if pp <= 0.15:
                return R, self.sa_dict[imc]['c3']
            elif pp < 0.8:
                c3 = self.sa_dict[imc]['c3']
                c4 = self.sa_dict[imc]['c4']
                return R, c3 + (c4 - c3) * log(pp / 0.15) / log(0.8 / 0.15)
            elif pp <= 5.0:
                return R, self.sa_dict[imc]['c4']
            else:
                # Not sure what's right here; should probably raise an error
                # but for now let's just use c4
                return R, self.sa_dict[imc]['c4']
        else:
            raise ValueError('unknown IMT %r' % imt)
    
