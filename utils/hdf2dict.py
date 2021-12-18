#! /usr/bin/env python

import h5py
import numpy as np

#
# You can use this file in two ways: You can call it as a program from
# a directory with a "shake_result.hdf" file in it, and it will spit out
# some slightly interesting information. You might do that as a test. Or
# you can put it in a directory with your code and do:
#
# from hdf2dict import get_result_dict
#
# and then call get_result_dict() with the path to a shake_result.hdf file
# as an argument, and it will return a dictionary with stuff in it.
#


def get_result_dict(filename):
    """
    This function is called with the path to a shake_result.hdf file. It
    returns a dictionary with the following keys:
        'mean': A grid of the conditional mean MMI
        'mean_metadata': A metadata dictionary for the mean grid
        'std': A grid of the conditional total standard deviation of MMI
        'std_metadata': A metadata dictionary for the std grid
        'phi': A grid of the prior within-event standard deviation of MMI
        'phi_metadata': A metadata dictionary for the phi grid
        'tau': A grid of the conditional between-event standard deviation of MMI
        'tau_metadata': A metadata dictionary for the tau grid
    The metadata dictionaries should be all the same. They have the following
    keys:
        'digits': some unknown thing
        'dx': the grid spacing in the X dimension
        'dy': the grid spacing in the Y dimension
        'nx': the number of grid points in the X dimension
        'ny': the number of grid points in the Y dimension
        'units': the physical units of the grids (intensity)
        'xmax': the maximum longitude of the grid
        'xmin': the minimum longitude of the grid
        'ymax': the maximum latitude of the grid
        'ymin': the minimum latitude of the grid
    The grids themselves will be of dimension (ny x nx) (rows x cols).
    """
    ddict = {}
    hdfobj = h5py.File(filename, "r+")
    group = hdfobj["arrays"]["imts"]["GREATER_OF_TWO_HORIZONTAL"]["MMI"]
    for name in ("mean", "std", "phi", "tau"):
        dset = group[name]
        ddict[name] = dset[()].astype(np.float64)
        metadata = {}
        for key, value in dset.attrs.items():
            metadata[key] = value
        ddict[name + "_metadata"] = metadata

    hdfobj.close()
    return ddict


if __name__ == "__main__":

    ddict = get_result_dict("shake_result.hdf")

    print("ddict keys: ", ddict.keys())
    print()
    print("Metadata:")
    print(ddict["mean_metadata"].items())
    print()
    print("Array shape: ", ddict["mean"].shape)
