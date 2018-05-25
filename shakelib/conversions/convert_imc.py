# Standard library imports
from abc import ABC, abstractmethod


class ComponentConverter(ABC):
    """Base class for implementing conversions between components."""

    def checkUnknown(self):
        """
        Checks if imc_in or imc_out is an unknown type.

        If of an unknown type, it is automatically set
        to AVERAGE_HORIZONTAL.
        """
        if self.imc_in not in self.conversion_graph:
            self.imc_in = 'Average horizontal'
        if self.imc_out not in self.conversion_graph:
            self.imc_out = 'Average horizontal'

    def convertAmps(self, imt, amps, rrups=None, mag=None):
        """
        Return an array of amps converted from one IMC to another.

        Note: This can be used to perform chained conversions.

        Args:
            imt (OpenQuake IMT): The intensity measure type of the input
                ground motions. Valid IMTs are PGA, PGV, and SA.
            amps (array): A numpy array of the (logged) ground motions
                to be converted.
            rrups (array): A numpy array of the same shape as amps,
                containing the rupture distances of the ground motions.
                Default is None.
            mag (float): The earthquake magnitude. Default is None.

        Returns:
            array: A numpy array of converted ground motions.
        """
        for idx, path in enumerate(self.path):
            if idx != len(self.path) - 1:
                self.imc_in = path
                self.imc_out = self.path[idx + 1]
                amps = self.convertAmpsOnce(imt, amps, rrups, mag)
        return amps

    @abstractmethod
    def convertAmpsOnce(self, imt, amps, rrups=None, mag=None):
        """
        Return an array of amps converted from one IMC to another.

        Note: This does not implement chained conversions.

        Args:
            imt (OpenQuake IMT): The intensity measure type of the input
                ground motions. Valid IMTs are PGA, PGV, and SA.
            amps (array): A numpy array of the (logged) ground motions
                to be converted.
            rrups (array): A numpy array of the same shape as amps,
                containing the rupture distances of the ground motions.
            mag (float): The earthquake magnitude.

        Returns:
            array: A numpy array of converted ground motions.
        """
        pass

    def convertSigmas(self, imt, sigmas):
        """
        Return an array of standard deviations converted from one IMC
        to another.

        Note:  This can be used to perform chained conversions.

        Args:
            imt (OpenQuake IMT): The intensity measure type of the input
                ground motions. Valid IMTs are PGA, PGV, and SA.
            sigmas (array): A numpy array of the standard deviations of
                the logged ground motions.

        Returns:
            array: A numpy array of converted standard deviations.
        """
        for idx, path in enumerate(self.path):
            if idx != len(self.path) - 1:
                self.imc_in = path
                self.imc_out = self.path[idx + 1]
                sigmas = self.convertSigmasOnce(imt, sigmas)
        return sigmas

    @abstractmethod
    def convertSigmasOnce(self, imt, sigmas):
        """
        Return an array of standard deviations converted from one IMC
        to another.

        Note: This does not implement chained conversions.

        Args:
            imt (OpenQuake IMT): The intensity measure type of the input
                ground motions. Valid IMTs are PGA, PGV, and SA.
            sigmas (array): A numpy array of the standard deviations of
                the logged ground motions.

        Returns:
            array: A numpy array of converted standard deviations.
        """
        pass

    @staticmethod
    def pathSearch(graph, imc_in, imc_out):
        """
        Helper method to create a "path" to convert one IMC to another.

        This can be used for chain conversions

        Args:
            graph (Dictionary): Dictionary of sets describing possible
                conversions between IMCs.
            imc_in (IMC): OpenQuake IMC type of the input amp array.
            imc_out (IMC): Desired OpenQuake IMC type of the output amps.

        Returns:
            list: IMCs as a path for to convert one IMC to another.
        """
        queue = [(imc_in, [imc_in])]
        while queue:
            (key, path) = queue.pop(0)
            for next_path in graph[key] - set(path):
                if next_path == imc_out:
                    yield path + [next_path]
                else:
                    queue.append((next_path, path + [next_path]))

    def getShortestPath(self, graph, imc_in, imc_out):
        """
        Create a "path" to convert one IMC to another.

        This can be used for chain conversions

        Args:
            graph (Dictionary): Dictionary of sets describing possible
                conversions between IMCs.
            imc_in (IMC): OpenQuake IMC type of the input amp array.
            imc_out (IMC): Desired OpenQuake IMC type of the output amps.

        Returns:
            list: IMCs as a path for to convert one IMC to another.

        Raises:
            ValueError if no path is found.
        """
        try:
            return next(self.pathSearch(graph, imc_in, imc_out))
        except StopIteration:
            if imc_in == imc_out:
                return [imc_in, imc_out]
            else:
                raise ValueError('No possible conversion between %r and %r.' %
                (imc_in, imc_out))

    @abstractmethod
    def _verifyConversion(self, imc_in, imc_out=None):
        """
        Helper method to ensure that the conversion is possible.

        Args:
            imc_in (IMC): OpenQuake IMC type of the input amp array.
            imc_out (IMC): Desired OpenQuake IMC type of the output amps.
                Default is None.

        Raises:
            ValueError if imc_in or imc_out are not valid..
        """
        pass
