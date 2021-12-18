#! /usr/bin/env python

import json
import numpy as np

#
# You can use this file in two ways: You can call it as a program from
# a directory with a "stationlist.json" file in it, and it will spit out
# some slightly interesting information. You might do that as a test. Or
# you can put it in a directory with your code and do:
#
# from json2dict import get_station_dict
#
# and then call get_station_dict() with the path to a stationlist.json file
# as an argument, and it will return a dictionary with stuff in it.
#


def get_station_dict(filename):
    """
    This function is called with the path to a stationlist.json file. It
    returns a dictionary with the following keys:
        'lons': An array of station longitudes
        'lats': An array of station latitudes
        'intensity': An array of observed intensities
        'residuals': An array of residuals (observed - predicted) (not
                     including bias)
        'phi': An array of within-event stddev of the predictions
        'tau': An array of between-event stddev of the predictions
        'extra_sigma': An array of the additional stddev of the observations
        'bias': An array of the bias of the data at the observation points
        'bias_stddev': An array of the stddev of the bias at the observation
                       points
    """
    lonlist = []
    latlist = []
    intensitylist = []
    residuallist = []
    philist = []
    taulist = []
    extra_sigma_list = []
    bias_list = []
    bias_sigma_list = []
    with open(filename) as f:
        sdata = json.load(f)

    for station in sdata["features"]:
        properties = station["properties"]
        if not properties["source"] == "DYFI":
            continue
        if not properties["intensity_flag"] == "0":
            continue
        intensity = properties["intensity"]
        if intensity is None or intensity == "null":
            continue
        intensitylist.append(float(intensity))
        extra_sigma_list.append(float(properties["intensity_stddev"]))
        for pred in properties["predictions"]:
            if not pred["name"] == "mmi":
                continue
            residuallist.append(float(intensity) - float(pred["value"]))
            philist.append(float(pred["phi"]))
            taulist.append(float(pred["tau"]))
            bias_list.append(float(pred["bias"]))
            bias_sigma_list.append(float(pred["bias_sigma"]))
            break
        lonlist.append(float(station["geometry"]["coordinates"][0]))
        latlist.append(float(station["geometry"]["coordinates"][1]))

    sdict = {
        "lons": np.array(lonlist),
        "lats": np.array(latlist),
        "intensity": np.array(intensitylist),
        "residuals": np.array(residuallist),
        "phi": np.array(philist),
        "tau": np.array(taulist),
        "extra_sigma": np.array(extra_sigma_list),
        "bias": np.array(bias_list),
        "bias_sigma": np.array(bias_sigma_list),
    }

    return sdict


if __name__ == "__main__":

    sdict = get_station_dict("stationlist.json")

    print("sdict keys: ", sdict.keys())
    print()
    print("Number of stations:")
    nsta = len(sdict["lons"])
    print(nsta)
    print()
    print("Min lon, Max lon, Min lat, Max lat:")
    if nsta > 0:
        print(
            np.min(sdict["lons"]),
            np.max(sdict["lons"]),
            np.min(sdict["lats"]),
            np.max(sdict["lats"]),
        )
    else:
        print("No data")
