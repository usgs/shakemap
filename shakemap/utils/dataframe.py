# stdlib imports
import json
import re
import string
import time

# third party imports
import pandas as pd
import numpy as np
from lxml import etree
from openpyxl import load_workbook, utils

COMPONENTS = ["GREATER_OF_TWO_HORIZONTALS", "GEOMETRIC_MEAN", "ARITHMETIC MEAN"]
CHANNEL_PATTERNS = [
    "^[H,B][H,L,N][E,N,Z,1,2,3]$",  # match standard seed names
    "^H[1,2]$",  # match H1/H2
    "^Z$",
]  # match Z
PGM_COLS = ["PGA", "PGV", "SA(0.3)", "SA(1.0)", "SA(3.0)"]
OPTIONAL = [
    "NAME",
    "DISTANCE",
    "REFERENCE",
    "INTENSITY",
    "SOURCE",
    "LOC",
    "INSTTYPE",
    "ELEV",
    "NRESP",
    "STDDEV",
    "",
    "FLAG",
    "INSTRUMENT",
    "PERIOD",
    "SENSITIVITY",
    "SERIAL",
    "SOURCE_FORMAT",
    "STRUCTURE",
    "DAMPING",
]
FLOATRE = "[-+]?[0-9]*\.?[0-9]+"


def dataframe_to_xml(df, xmlfile, reference=None):
    """Write an MMI dataframe to ShakeMap XML format.

    This method accepts a dataframe with this structure:
     - STATION: Station code (REQUIRED)
     - LAT: Station latitude. (REQUIRED)
     - LON: Station longitude. (REQUIRED)
     - DISTANCE: Distance (km) from station to origin.
     - NETID: Network ID
     - SOURCE: Description of data contributor.
     - INTENSITY: MMI intensity.
     - NRESP: Number of responses for aggregated intensity.

    Args:
        df (DataFrame): Pandas dataframe, as described in read_excel.
        xmlfile (str): Path to file where XML file should be written.
    """
    root = etree.Element("shakemap-data", code_version="3.5", map_version="3")

    create_time = int(time.time())
    stationlist = etree.SubElement(root, "stationlist", created=f"{int(create_time):d}")
    if reference is not None:
        stationlist.attrib["reference"] = reference

    processed_stations = []

    for _, row in df.iterrows():
        tmprow = row.copy()

        # assign required columns
        stationcode = str(tmprow["STATION"]).strip()

        netid = tmprow["NETID"].strip()
        if not stationcode.startswith(netid):
            stationcode = f"{netid}.{stationcode}"

        # if this is a dataframe created by shakemap,
        # there will be multiple rows per station.
        # below we process all those rows at once,
        # so we need this bookkeeping to know that
        # we've already dealt with this station
        if stationcode in processed_stations:
            continue

        station = etree.SubElement(stationlist, "station")

        station.attrib["code"] = stationcode
        station.attrib["lat"] = f"{tmprow['LAT']:.4f}"
        station.attrib["lon"] = f"{tmprow['LON']:.4f}"

        # assign optional columns
        if "NETID" in tmprow:
            station.attrib["netid"] = tmprow["NETID"].strip()
        if "DISTANCE" in tmprow:
            station.attrib["dist"] = f"{tmprow['DISTANCE']:.1f}"
        if "INTENSITY" in tmprow:
            station.attrib["intensity"] = f"{tmprow['INTENSITY']:.1f}"
        if "STDDEV" in tmprow:
            station.attrib["intensity_stddev"] = f"{tmprow['STDDEV']:.4f}"
        if "NRESP" in tmprow:
            station.attrib["nresp"] = f"{int(tmprow['NRESP']):d}"
        if "SOURCE" in tmprow:
            station.attrib["source"] = tmprow["SOURCE"].strip()

        processed_stations.append(stationcode)

    tree = etree.ElementTree(root)
    tree.write(xmlfile, pretty_print=True)
