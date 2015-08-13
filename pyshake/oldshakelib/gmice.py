#!/usr/bin/env python

class GMICE(object):
    def __init__(self,magnitude):
        """
        Instantiate a GMICE object.
        @param magnitude: Earthquake magnitude.
        """
        pass

    def convertToIntensity(self,pga,pgv):
        """
        Convert PGA/PGV pairs to MMI using GMICE specific algorithm.
        Implemented in sub-class.
        @param pga: Numpy array of PGA values.
        @param pgv: Numpy array of PGV values.
        @return: Numpy array of MMI values.
        """
        pass

    def convertFromIntensity(self,mmi):
        """
        Convert MMI values to PGA/PGV.
        Implemented in sub-class.
        @param mmi: Numpy array of MMI values.
        @return: Tuple of two numpy arrays of (pga,pgv) values.
        """
        pass

    #static methods - that is, methods that are associated with the class but do not require an 
    #instantiated object.
    @staticmethod
    def uncertaintyFunction():
        """
        A static class method implemented in sub-classes, providing a function for calculating uncertainties.
        This will be used both by the sub-class itself, and the Data object for calculating uncertainties.
        NB: Not entirely sure about the syntax of the @staticmethod decorator... putting it in the base class
        may not work.  Experimentation required.
        """
        pass
