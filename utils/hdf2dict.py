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
    In addition to the above grids, the following arrays will be returned:
        'Sigma_HH_YD': A numpy array, see below for possible values.
        'Sigma_HH_YD_metadata': This will always be 'None'.
        'C': A numpy array, see below for possible values.
        'C_metadata': This will always be 'None'.
        'add_uncertainty': A numpy array, see below for possible values.
        'add_uncertainty_metadata': This will always be 'None'.
        'sta_per_ix': A numpy array, see below for possible values.
        'sta_per_ix_metadata': This will always be 'None'.
    These arrays will have the following contents:
        'None': If this function is run on an older version of shake_results.hdf
                that was produced by model.py before it generated the output
                needed to produce these arrays.
        '[]': I.e., empty arrays. When the shakemap contained no input station data
              sufficient to generate the arrays.
        ...: Arrays of shapes appropriate to the input and output data dimensions.
             It is incumbent upon the user to make proper use of these arrays.
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
    for name in (
        "add_uncertainty",
        "Sigma_HH_YD",
        "C",
        "sta_per_ix",
        "sta_phi",
        "sta_lons_rad",
        "sta_lats_rad",
    ):
        ddict[name + "_metadata"] = None
        try:
            dset = group[name]
            ddict[name] = dset[()].astype(np.float64)
        except Exception:
            ddict[name] = None

    hdfobj.close()
    return ddict


if __name__ == "__main__":

    ddict = get_result_dict("shake_result.hdf")

    print("ddict keys: ", ddict.keys())
    print()
    for key in ddict.keys():
        if ddict[key] is None:
            print(f"ddict {key} is None")
    for name in ("mean", "std", "phi", "tau"):
        print(f"{name} shape is {(ddict[name].shape)}")
    for key in (
        "add_uncertainty",
        "Sigma_HH_YD",
        "C",
        "sta_per_ix",
        "sta_phi",
        "sta_lons_rad",
        "sta_lats_rad",
    ):
        if ddict[key] is None:
            print(f"{key} shape is None")
        else:
            print(f"{key} shape is {(ddict[key].shape)}")
    print("Metadata:")
    print(ddict["mean_metadata"].items())
    print()
    print("Array shape: ", ddict["mean"].shape)
