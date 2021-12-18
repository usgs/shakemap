#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
from collections import OrderedDict
import pandas as pd

from libcomcat.search import get_event_by_id


def main():
    desc = "Program to get CSV file of station data from ShakeMap."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("eventid", help="Comcat event ID.")
    args = parser.parse_args()
    evid = args.eventid

    # Download stationlist
    event = get_event_by_id(evid)
    shakemap = event.getProducts("shakemap")[0]
    json_file = evid + "_stationlist.json"
    shakemap.getContent("stationlist.json", json_file)
    with open(json_file) as f:
        station_dict = json.load(f)

    # Extract info in tabular form
    out_dict = OrderedDict()
    out_dict["lat"] = []
    out_dict["lon"] = []
    out_dict["rjb"] = []
    out_dict["repi"] = []
    out_dict["pga_percent_g"] = []
    out_dict["pgv_cm_s"] = []
    for f in station_dict["features"]:
        if f["properties"]["station_type"] == "seismic":
            out_dict["lon"].append(f["geometry"]["coordinates"][0])
            out_dict["lat"].append(f["geometry"]["coordinates"][1])
            out_dict["rjb"].append(f["properties"]["distances"]["rjb"])
            out_dict["repi"].append(f["properties"]["distances"]["repi"])
            out_dict["pga_percent_g"].append(f["properties"]["pga"])
            out_dict["pgv_cm_s"].append(f["properties"]["pgv"])

    out_file = evid + "_stationlist.csv"
    out_df = pd.DataFrame(out_dict)
    out_df.to_csv(out_file, index=False)


if __name__ == "__main__":
    main()
