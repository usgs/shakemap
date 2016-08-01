
class BommerAlarcon2006(object):

    """
    Class for conversion between PGV (units of cm/s) and PSA05 (units of g)
    by Alarcon and Bommer (2006).

    Comments for testing.
    """

    __vfact = 1.0 / (20.0) * 100.0 * 9.81

    @classmethod
    def pgv2psa05(cls, pgv):
        """
        Convert PGV in cm/s to PSA05 in g.
        ** PGV must be linear units **
        :param pgv:
            Numpy array or float of PGV values; linear units.
        : returns:
            Values converted to PSA05.
        """
        return pgv / cls.__vfact

    @classmethod
    def psa052pgv(cls, psa05):
        """
        Convert PSA05 in g to PGV cm/s.
        ** PSA10 must be linear units. **
        :param psa05:
            Numpy array or float of PSA05 values; linear units.
        : returns:
            Values converted to PGV.
        """
        return psa05 * cls.__vfact

    def getVfact(self):
        return self.__vfact
