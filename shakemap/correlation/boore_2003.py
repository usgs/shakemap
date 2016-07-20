import numpy as np


class Boore2003(object):
    def getSpatialCorrelation(self, dists, imt):
        return 1.0 - np.exp(-1.0 * np.sqrt(0.6 * dists))

