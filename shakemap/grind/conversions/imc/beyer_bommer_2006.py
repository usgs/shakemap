
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, PGV, SA

import numpy as np

class BeyerBommer2006(object):
    """
    Implements conversion for various "Intensity Measure
    Components" (IMCs) per Beyer and Bommer (2006).
    
    IMC equivalencies:
    
    +---------------------------+------------------------+
    | OpenQuake                 | Beyer & Bommer         |
    +===========================+========================+
    | AVERAGE_HORIZONTAL        | GMxy (Geometric mean)  |
    +---------------------------+------------------------+
    | HORIZONTAL                | ???                    |
    +---------------------------+------------------------+
    | MEDIAN_HORIZONTAL         | AMxy (Arithmetic mean) |
    +---------------------------+------------------------+
    | GMRotI50                  | GMRotI50               |
    +---------------------------+------------------------+
    | RotD50                    | GMRotD50               |
    +---------------------------+------------------------+
    | RANDOM_HORIZONTAL         | Random                 |
    +---------------------------+------------------------+
    | GREATER_OF_TWO_HORIZONTAL | Env_xy                 |
    +---------------------------+------------------------+
    | VERTICAL                  | ---                    |
    +---------------------------+------------------------+

    Note that AVERAGE_HORIZONTAL is the "reference" type

    To do
        - Inherit from ConvertIMC class. 
    
    References: 
        Beyer, K., & Bommer, J. J. (2006). Relationships between median values
        and between aleatory variabilities for different definitions of the 
        horizontal component of motion. Bulletin of the Seismological Society of
        America, 96(4A), 1512-1522. 
        `[link] <http://www.bssaonline.org/content/96/4A/1512.short>`__
    """
    __pga_pgv_col_names = ['c12', 'c34', 'R']
    __sa_col_names = ['c1', 'c2', 'c3', 'c4', 'R']
    __pga_dict = {
        const.IMC.MEDIAN_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.01, 1.00]))),
        const.IMC.GMRotI50: dict(list(zip(__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
        const.IMC.RotD50: dict(list(zip(__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
        const.IMC.RANDOM_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
        const.IMC.GREATER_OF_TWO_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.1, 0.05, 1.02])))}
    __pgv_dict = {
        const.IMC.MEDIAN_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.01, 1.00]))),
        const.IMC.GMRotI50: dict(list(zip(__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
        const.IMC.RotD50: dict(list(zip(__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
        const.IMC.RANDOM_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
        const.IMC.GREATER_OF_TWO_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.1, 0.05, 1.02])))}
    __sa_dict = {
        const.IMC.MEDIAN_HORIZONTAL: dict(list(zip(__sa_col_names, [1.0, 1.0, 0.01, 0.02, 1.00]))),
        const.IMC.GMRotI50: dict(list(zip(__sa_col_names, [1.0, 1.0, 0.03, 0.04, 1.00]))),
        const.IMC.RotD50: dict(list(zip(__sa_col_names, [1.0, 1.0, 0.02, 0.03, 1.00]))),
        const.IMC.RANDOM_HORIZONTAL: dict(list(zip(__sa_col_names, [1.0, 1.0, 0.07, 0.11, 1.05]))),
        const.IMC.GREATER_OF_TWO_HORIZONTAL: dict(list(zip(__sa_col_names, [1.1, 1.2, 0.04, 0.07, 1.02])))}
    
    @staticmethod
    def ampIMCtoIMC(amps, imc_in, imc_out, imt):
        """ 
        Returns amps converted from one IMC to another.

        **Important**: 

            - Assumes the input amps are in linear (not log) space 
            - IMC types 'VERTICAL' and 'HORIZONTAL' are not supported
        
        :param amps: 
            Numpy array of ground motion amplitudes. 
        :param imc_in:
            OpenQuake IMC type of the input amp array. 
            `[link] <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__
        :param imc_out:
            Desired OpenQuake IMC type of the output amps. 
            `[link] <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__
        :param imt:
            OpenQuake IMT of the input amps (must be one of PGA, PGV, or SA). 
            `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`
        :returns:
            Numpy array of amps converted from imc_in to imc_out
        """
        if imc_in == const.IMC.AVERAGE_HORIZONTAL:
            # The amps are already in the B&B "reference" type ("GM", i.e.,
            # geometric mean)
            denom = 1
        elif imc_in == const.IMC.GREATER_OF_TWO_HORIZONTAL or \
             imc_in == const.IMC.MEDIAN_HORIZONTAL or \
             imc_in == const.IMC.GMRotI50 or \
             imc_in == const.IMC.RotD50 or \
             imc_in == const.IMC.RANDOM_HORIZONTAL:
            denom = BeyerBommer2006.__GM2other(imt, imc_in)
        else:
            raise ValueError('unknown IMC %r' % imc_in)

        if imc_out == const.IMC.AVERAGE_HORIZONTAL:
            # The previous step will put the amps into the B&B "reference"
            # type ("GM", i.e. geometric mean)
            numer = 1
        elif imc_out == const.IMC.GREATER_OF_TWO_HORIZONTAL or \
             imc_out == const.IMC.MEDIAN_HORIZONTAL or \
             imc_out == const.IMC.GMRotI50 or \
             imc_out == const.IMC.RotD50 or \
             imc_out == const.IMC.RANDOM_HORIZONTAL:
            numer = BeyerBommer2006.__GM2other(imt, imc_out)
        else:
            raise ValueError('unknown IMC %r' % imc_out)
        
        return amps * (numer / denom)

    @staticmethod
    def sigmaIMCtoIMC(sigmas, imc_in, imc_out, imt):
        """ 
        Returns standard deviations converted from one IMC to another.

        **Important**: 

            - Assumes the input sigmas are in log space
            - IMC types 'VERTICAL' and 'HORIZONTAL' are not supported
        
        :param sigmas:
            Numpy array of standard deviations. 
        :param imc_in:
            OpenQuake IMC type of the input sigmas array. 
            `[link] <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__
        :param imc_out:
            Desired OpenQuake IMC type of the output sigmas. 
            `[link] <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__
        :param imt:
            OpenQuake IMT of the input sigmas (must be one of PGA, PGV, or SA)
            `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__

        :returns:
            Numpy array of standard deviations converted from imc_in to imc_out
        """
        if imc_in == const.IMC.AVERAGE_HORIZONTAL:
            # The amps are already in the B&B "reference" type ("GM", i.e.,
            # geometric mean)
            R = 1
            sig_log_ratio = 0
        elif imc_in == const.IMC.GREATER_OF_TWO_HORIZONTAL or \
             imc_in == const.IMC.MEDIAN_HORIZONTAL or \
             imc_in == const.IMC.GMRotI50 or\
             imc_in == const.IMC.RotD50 or\
             imc_in == const.IMC.RANDOM_HORIZONTAL:
            R, sig_log_ratio = BeyerBommer2006.__GM2otherSigma(imt, imc_in)
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
             imc_out == const.IMC.GMRotI50 or\
             imc_out == const.IMC.RotD50 or\
             imc_out == const.IMC.RANDOM_HORIZONTAL:
            R, sig_log_ratio = BeyerBommer2006.__GM2otherSigma(imt, imc_out)
        else:
            raise ValueError('unknown IMC %r' % imc_out)

        sigma_out = np.sqrt(sigma_GM2 * R**2 + sig_log_ratio**2)
        return sigma_out

    @staticmethod
    def __GM2other(imt, imc):
        """ Helper function to extract coefficients from the parameter tables """
        if 'PGA' in imt:
            return BeyerBommer2006.__pga_dict[imc]['c12']
        elif 'PGV' in imt:
            return BeyerBommer2006.__pgv_dict[imc]['c12']
        elif 'SA' in imt:
            pp = imt.period
            if pp <= 0.15:
                return BeyerBommer2006.__sa_dict[imc]['c1']
            elif pp < 0.8:
                c1 = BeyerBommer2006.__sa_dict[imc]['c1']
                c2 = BeyerBommer2006.__sa_dict[imc]['c2']
                return c1 + (c2 - c1) * np.log(pp / 0.15) / np.log(0.8 / 0.15)
            elif pp <= 5.0:
                return BeyerBommer2006.__sa_dict[imc]['c2']
            else:
                # Not sure what's right here; should probably raise an error
                # but for now let's just use c2
                return BeyerBommer2006.__sa_dict[imc]['c2']
        else:
            raise ValueError('unknown IMT %r' % imt)

    @staticmethod
    def __GM2otherSigma(imt, imc):
        """ Helper function to extract coefficients from the parameter tables """
        if 'PGA' in imt:
            return BeyerBommer2006.__pga_dict[imc]['R'], BeyerBommer2006.__pga_dict[imc]['c34']
        elif 'PGV' in imt:
            return BeyerBommer2006.__pgv_dict[imc]['R'], BeyerBommer2006.__pgv_dict[imc]['c34']
        elif 'SA' in imt:
            R = BeyerBommer2006.__sa_dict[imc]['R']
            pp = imt.period
            if pp <= 0.15:
                return R, BeyerBommer2006.__sa_dict[imc]['c3']
            elif pp < 0.8:
                c3 = BeyerBommer2006.__sa_dict[imc]['c3']
                c4 = BeyerBommer2006.__sa_dict[imc]['c4']
                return R, c3 + (c4 - c3) * np.log(pp / 0.15) / np.log(0.8 / 0.15)
            elif pp <= 5.0:
                return R, BeyerBommer2006.__sa_dict[imc]['c4']
            else:
                # Not sure what's right here; should probably raise an error
                # but for now let's just use c4
                return R, BeyerBommer2006.__sa_dict[imc]['c4']
        else:
            raise ValueError('unknown IMT %r' % imt)
