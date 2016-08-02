
class NewmarkHall1982(object):

    """
    Class for conversion between PGA and PSA10 by Newmark and Hall (1982)
    """

    __vfact = 37.27 * 2.54

    """
    37.27 is derived from ->   SA(f) = V * 2 * pi * f * 1.65 / 386.09
    where:
    SA(f) is spectral acceleration at frequency f (in g)
    f is the frequency of interest
    V is the velocity
    1.65 is the N&H amplification factor for velocity at 5% damping
    386.09 is the acceleration of gravity in inches per second per g

    2.54 is the conversion factor to convert from cm/s to in/s
    """

    @staticmethod
    def pgv2psa10(pgv):
        """
        Convert PGV in cm/s to PSA10 in g.
        ** PGV must be linear units **
        :param pgv:
            Numpy array or float of PGV values; linear units.
        : returns:
            Values converted to PSA10.
        """
        return pgv / NewmarkHall1982.__vfact

    @staticmethod
    def psa102pgv(psa10):
        """
        Convert PSA10 in g to PGV cm/s.
        ** PSA10 must be linear units **
        :param psa10:
            Numpy array or float of PSA10 values; linear units.
        : returns:
            Values converted to PGV.
        """
        return psa10 * NewmarkHall1982.__vfact

    def getVfact(self):
        return self.__vfact

