from abc import ABC, abstractmethod


class CrossCorrelationBase(ABC):
    """
    The abstract base class of cross-correlation functions.
    """

    @abstractmethod
    def __init__(self, periods):
        """The derived class should take an array of periods as its single
        argument. These are the periods to which the index arrays of the
        getCorrelation() method will refer.
        """
        pass  # pragma: no cover

    @abstractmethod
    def getCorrelation(self, ix1, ix2, h):
        """
        Compute the correlation between two periods and a separation distance
        of h km. ix1 and ix2 give the indices into the array of periods that
        is provided to the class constructor. The code is implemented in this
        somewhat awkward way for performance reasons. See the loth_baker_2013
        module for an example.

        The index arrays (ix1 and ix2) and h array must have the same
        dimensions. The
        indices may be equal, and there is no restriction on which one is
        larger. The indices refer to periods in the 'period' argument to the
        class constructor.

        Args:
            ix1 (ndarray):
                The indices of the first period of interest.
            ix2 (ndarrays):
                The indices of the second period of interest.
            h (ndarray):
                The separation distance between two sites (units of km).
                The result is placed in h, so h must be copied if it is
                to be preserved.

        Returns:
            h (ndarray): The predicted correlation coefficient. The output
            array will have the same shape as the inputs.

        """
        pass  # pragma: no cover
