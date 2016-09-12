
from openquake.hazardlib import const

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
    | HORIZONTAL                | Random?                |
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

    Notes 

        - AVERAGE_HORIZONTAL is the "reference" type. 
        - The OQ IMC "HORIZONAL" indicates that the horizontal IMC category
          may be ambiguous. In these cases, we are assuming that it is a 
          random horizontal component as a default.

    To do

        - Inherit from ConvertIMC class. 

    References

        Beyer, K., & Bommer, J. J. (2006). Relationships between median values
        and between aleatory variabilities for different definitions of the 
        horizontal component of motion. Bulletin of the Seismological Society of
        America, 96(4A), 1512-1522. 
        `[link] <http://www.bssaonline.org/content/96/4A/1512.short>`__
    """

    # c12 = median ratio
    # c34 = std. of log ratio
    # R = ratio of sigma_pga values
    __pga_pgv_col_names = ['c12', 'c34', 'R']
    __sa_col_names = ['c1', 'c2', 'c3', 'c4', 'R']
    __pga_dict = {
        const.IMC.MEDIAN_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.01, 1.00]))),
        const.IMC.GMRotI50: dict(list(zip(__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
        const.IMC.RotD50: dict(list(zip(__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
        const.IMC.RANDOM_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
        const.IMC.HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
        const.IMC.GREATER_OF_TWO_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.1, 0.05, 1.02])))}
    __pgv_dict = {
        const.IMC.MEDIAN_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.01, 1.00]))),
        const.IMC.GMRotI50: dict(list(zip(__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
        const.IMC.RotD50: dict(list(zip(__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
        const.IMC.RANDOM_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
        const.IMC.HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
        const.IMC.GREATER_OF_TWO_HORIZONTAL: dict(list(zip(__pga_pgv_col_names, [1.1, 0.05, 1.02])))}
    __sa_dict = {
        const.IMC.MEDIAN_HORIZONTAL: dict(list(zip(__sa_col_names, [1.0, 1.0, 0.01, 0.02, 1.00]))),
        const.IMC.GMRotI50: dict(list(zip(__sa_col_names, [1.0, 1.0, 0.03, 0.04, 1.00]))),
        const.IMC.RotD50: dict(list(zip(__sa_col_names, [1.0, 1.0, 0.02, 0.03, 1.00]))),
        const.IMC.RANDOM_HORIZONTAL: dict(list(zip(__sa_col_names, [1.0, 1.0, 0.07, 0.11, 1.05]))),
        const.IMC.HORIZONTAL: dict(list(zip(__sa_col_names, [1.0, 1.0, 0.07, 0.11, 1.05]))),
        const.IMC.GREATER_OF_TWO_HORIZONTAL: dict(list(zip(__sa_col_names, [1.1, 1.2, 0.04, 0.07, 1.02])))}

    @staticmethod
    def ampIMCtoIMC(amps, imc_in, imc_out, imt):
        """ 
        Returns amps converted from one IMC to another.

        **Important**:

            - Assumes the input amps are in natural log (not linear) space
            - IMC type 'VERTICAL' is not supported

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
            Numpy array of amps converted from imc_in to imc_out.
        """
        
        denom = BeyerBommer2006.__GM2other(imt, imc_in)
        numer = BeyerBommer2006.__GM2other(imt, imc_out)

        return amps + np.log(numer / denom)

    @staticmethod
    def sigmaIMCtoIMC(sigmas, imc_in, imc_out, imt):
        """ 
        Returns standard deviations converted from one IMC to another.

        **Important**: 

            - Assumes the input sigmas are in natural log space
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

        #---------------------------------------------------
        # Take input sigma to geometric mean sigma. 
        # Solve eqn 8 for sigma_GM2, which is sigma_logSa_GM
        # This is the sigma converted to geometric mean
        # (i.e., reference component). 
        #---------------------------------------------------
        R, sig_log_ratio = BeyerBommer2006.__GM2otherSigma(imt, imc_in)
        sigma_GM2 = (sigmas**2 - sig_log_ratio**2) / R**2

        #---------------------------------------------------
        # Evaluate equation 8 to go from GM to requested IMC
        #---------------------------------------------------
        R, sig_log_ratio = BeyerBommer2006.__GM2otherSigma(imt, imc_out)
        sigma_out = np.sqrt(sigma_GM2 * R**2 + sig_log_ratio**2)

        return sigma_out

    @staticmethod
    def __checkIMC(imc):
        """
        Check that the requested IMC is available.
        """
        if imc != const.IMC.GREATER_OF_TWO_HORIZONTAL and \
           imc != const.IMC.MEDIAN_HORIZONTAL and \
           imc != const.IMC.GMRotI50 and \
           imc != const.IMC.RotD50 and \
           imc != const.IMC.RANDOM_HORIZONTAL and \
           imc != const.IMC.HORIZONTAL:
            raise ValueError('unknown IMC %r' % imc)


    @staticmethod
    def __GM2other(imt, imc):
        """
        Helper function to extract coefficients from the parameters for 
        converting the median ground motions.

        :param imt:
            Intensity measure type.
        :param imc:
            Intensity measure component.
        :returns:
            Median ratios. This is directly from Table 2 for PGA and PGV, and
            computed from coefficients in Table 3 along with eqn 10 for Sa. 
        """

        if imc == const.IMC.AVERAGE_HORIZONTAL:
            # The amps are already in the B&B "reference" type ("GM", i.e.,
            # geometric mean)
            return 1.

        BeyerBommer2006.__checkIMC(imc)


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
        """
        Helper function to extract coefficients from the parameters for
        converting standard deviations. 

        :param imt:
            Intensity measure type.
        :param imc:
            Intensity measure component.
        :returns:
            Coefficients: R (ratio of sigma values), standard deviation of sigma ratios.  
        """

        if imc == const.IMC.AVERAGE_HORIZONTAL:
            # The amps are already in the B&B "reference" type ("GM", i.e.,
            # geometric mean)
            return 1., 0.

        BeyerBommer2006.__checkIMC(imc)

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
