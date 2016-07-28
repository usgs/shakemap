#!/usr/bin/env python

import numpy as np


class Zone(object):

    def __init__(self):
        self.zonedict = {}

    def addZone(self, name, depths, magnitudes, gmpes):
        ndep = len(depths)
        nmag = len(magnitudes)
        if not isinstance(gmpes, np.ndarray):
            raise TypeError(
                'Input gmpe list must be a numpy array of dimensions: %i %i' % (ndep, nmag))
        depths = np.array(depths)
        magnitudes = np.array(magnitudes)
        # make sure that depths and magnitudes are increasing
        if np.any(np.diff(depths) < 0) or np.any(np.diff(magnitudes) < 0):
            raise ValueError(
                'Depths and magnitudes must be monotonically increasing')
