import numpy as np

from openquake.hazardlib import const


class AbrahamsonBhasin2020(object):
    """
    Module to make conditional conversions from SA(T) to PGV. The user
    provides a magnitude when an instance and the instance can report back
    what period of SA is appropriate to use for computing PGV. The module
    is based upon:

    Abrahamson, N.A., and Bhasin, S. (2020) Conditional ground-motion model
    for peak ground velocity for active crustal regions. PEER Report 2020/05,
    Pacific Earthquake Engineering Research Center, U.C. Berkeley.
    """

    def __init__(self, mag):
        self.mag = mag
        mag_T = [
            (3.5, 0.20),
            (4.5, 0.28),
            (5.5, 0.40),
            (6.5, 0.95),
            (7.5, 1.40),
            (8.5, 2.80),
        ]
        self.Tref = 0
        for mm, tt in mag_T:
            if np.any(mag <= mm):
                self.Tref = tt
                break
        if self.Tref == 0:
            self.Tref = tt
        self.coeff = {
            "a1": 5.39,
            "a2": 0.799,
            "a3": 0.654,
            "a4": 0.479,
            "a5": -0.062,
            "a6": -0.359,
            "a7": -0.134,
            "a8": 0.023,
            "phi": 0.29,
            "tau": 0.16,
            "sigma": 0.33,
        }

    def getTref(self):
        """
        Returns the period of the SA amps suitable for conversion to PGV.
        Args:
            None.

        Returns:
            (float): The period of the SA amps to supply to getPGVandSTDDEVS()
        """
        return self.Tref

    def getPGVandSTDDEVS(self, psa, stddevs, sdtypes, rrup, vs30):
        """
        Function to compute the PGVs and inflate the uncertainties based on
        the input SA amps and uncertainties.

        Args:
            psa (np.array): An array os SA amp of the appropriate period.
            stddevs(np.arrays): A list of one or more array of standard
                deviations of the psa amps. The list may also contain arrays
                of the stddevs inflated by the point-source adjustments. Each
                array must be the same shape as the psa array.
            sdtypes (list): A list of the types of the standard deviations in
                stddevs (in order). These are OpenQuake const.StdDev types.
            rrup (np.array): Rupture distance to the points in psa. Must be
                the same shape as psa.
            vs30 (np.array): The Vs30 value at the points in psa. Must be the
                same shape as psa.

        Returns:
            (tuple: ndarray, list): An array of PGV values, and a list of
            stddev arrays corresponding to the input stddev types, now
            with the uncertainty propagated through the conversion.
        """
        c = self.coeff
        m = self.mag
        if np.any(m < 5):
            f1 = c["a2"]
        elif np.any(m <= 7.5):
            f1 = c["a2"] + (c["a3"] - c["a2"]) * (m - 5.0) / 2.5
        else:
            f1 = c["a3"]

        pgv = (
            c["a1"]
            + f1 * psa
            + c["a4"] * (m - 6.0)
            + c["a5"] * (8.5 - m) ** 2
            + c["a6"] * np.log(rrup + 5.0 * np.exp(0.4 * (m - 6.0)))
            + (c["a7"] + c["a8"] * (m - 5.0)) * np.log(vs30 / 425)
        )

        #
        # ShakeMap will often produce lists of stddevs that are twice
        # as long as sdtypes: one set that has been inflated by the
        # point-source to finite-fault factors, and one that has not.
        # Here we modify both, if they are present.
        #
        stddevs_out = []
        mytypes = sdtypes.copy()
        if len(stddevs) == 2 * len(sdtypes):
            mytypes.extend(sdtypes)
        for i, sdtype in enumerate(mytypes):
            if sdtype == const.StdDev.INTER_EVENT:
                sdout = np.sqrt(f1 ** 2 * stddevs[i] ** 2 + c["tau"] ** 2)
            elif sdtype == const.StdDev.INTRA_EVENT:
                sdout = np.sqrt(f1 ** 2 * stddevs[i] ** 2 + c["phi"] ** 2)
            else:
                sdout = np.sqrt(f1 ** 2 * stddevs[i] ** 2 + c["sigma"] ** 2)
            stddevs_out.append(sdout)

        return pgv, stddevs_out


class AbrahamsonBhasin2020PGA(AbrahamsonBhasin2020):
    """
    An alternative to AbrahamsonBhasin2020 when only PGA is available for
    conversion.
    """

    def __init__(self, mag):
        self.mag = mag
        self.Tref = 0
        self.coeff = {
            "a1": 4.77,
            "a2": 0.738,
            "a3": 0.484,
            "a4": 0.275,
            "a5": -0.036,
            "a6": -0.332,
            "a7": -0.44,
            "a8": 0.0,
            "phi1": 0.32,
            "phi2": 0.42,
            "phi": 0.0,
            "tau1": 0.12,
            "tau2": 0.26,
            "tau": 0.0,
            "sigma1": 0.34,
            "sigma2": 0.49,
            "sigma": 0.0,
            "m1": 5.0,
            "m2": 7.0,
        }
        c = self.coeff
        mref = (self.mag - c["m1"]) / (c["m2"] - c["m1"])

        #
        # Here we set the magnitude-dependent phi, tau, and sigma
        # to be used in the parent cless' getPGVandSTDDEVS()
        #
        if self.mag < c["m1"]:
            c["phi"] = c["phi1"]
            c["tau"] = c["tau1"]
            c["sigma"] = c["sigma1"]
        elif self.mag <= c["m2"]:
            c["phi"] = c["phi1"] + (c["phi2"] - c["phi2"]) * mref
            c["tau"] = c["tau1"] + (c["tau2"] - c["tau2"]) * mref
            c["sigma"] = c["sigma1"] + (c["sigma2"] - c["sigma2"]) * mref
        else:
            c["phi"] = c["phi2"]
            c["tau"] = c["tau2"]
            c["sigma"] = c["sigma2"]


class AbrahamsonBhasin2020SA1(AbrahamsonBhasin2020):
    """
    An alternative to AbrahamsonBhasin2020 when only SA(1.0) is available for
    conversion.
    """

    def __init__(self, mag):
        self.mag = mag
        self.Tref = 1.0
        self.coeff = {
            "a1": 4.80,
            "a2": 0.82,
            "a3": 0.55,
            "a4": 0.27,
            "a5": 0.054,
            "a6": -0.382,
            "a7": -0.21,
            "a8": 0.0,
            "phi1": 0.28,
            "phi2": 0.38,
            "phi": 0.0,
            "tau1": 0.12,
            "tau2": 0.17,
            "tau": 0.0,
            "sigma1": 0.30,
            "sigma2": 0.42,
            "sigma": 0.0,
            "m1": 5.0,
            "m2": 7.0,
        }
        c = self.coeff
        mref = (self.mag - c["m1"]) / (c["m2"] - c["m1"])

        #
        # Here we set the magnitude-dependent phi, tau, and sigma
        # to be used in the parent cless' getPGVandSTDDEVS()
        #
        if self.mag < c["m1"]:
            c["phi"] = c["phi1"]
            c["tau"] = c["tau1"]
            c["sigma"] = c["sigma1"]
        elif self.mag <= c["m2"]:
            c["phi"] = c["phi1"] + (c["phi2"] - c["phi2"]) * mref
            c["tau"] = c["tau1"] + (c["tau2"] - c["tau2"]) * mref
            c["sigma"] = c["sigma1"] + (c["sigma2"] - c["sigma2"]) * mref
        else:
            c["phi"] = c["phi2"]
            c["tau"] = c["tau2"]
            c["sigma"] = c["sigma2"]
