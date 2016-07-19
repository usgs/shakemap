
import pandas as pd


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
    
    @classmethod
    def pgv2psa10(cls,pgv):
        """
        Convert PGV in cm/s to PSA10 in g.
        ** PGV must be linear units **
        :param pgv:
            Numpy array or float of PGV values; linear units.
        : returns:
            Values converted to PSA10. 
        """
        return pgv/cls.__vfact

    @classmethod
    def psa102pgv(cls,psa10):
        """
        Convert PSA10 in g to PGV cm/s.
        ** PSA10 must be linear units **
        :param psa10:
            Numpy array or float of PSA10 values; linear units.
        : returns:
            Values converted to PGV. 
        """
        return psa10*cls.__vfact


class BommerAlarcon2006(object):
    
    """
    Class for conversion between PGA and PSA10 by Beyer and Bommer (2006).

    Comments for testing.
    """
    
    __vfact = 20.0 * 9.81

    @classmethod
    def pgv2psa10(cls,pgv):
        """
        Convert PGV in cm/s to PSA10 in g.
        ** PGV must be linear units **
        :param pgv:
            Numpy array or float of PGV values; linear units.
        : returns:
            Values converted to PSA10. 
        """
        return pgv/cls.__vfact

    @classmethod
    def psa102pgv(cls,psa10):
        """
        Convert PSA10 in g to PGV cm/s.
        ** PSA10 must be linear units. **
        :param psa10:
            Numpy array or float of PSA10 values; linear units.
        : returns:
            Values converted to PGV. 
        """
        return psa10*cls.__vfact
