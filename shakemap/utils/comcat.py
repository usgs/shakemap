# stdlib imports
import json
import os.path
from io import StringIO
from urllib.request import urlopen

DETAIL_TEMPLATE = (
    "https://earthquake.usgs.gov/fdsnws/event/1/query?eventid=[EVENTID]&format=geojson"
)


def get_detail_json(eventid):
    """
    Return the detailed JSON dictionary for a ComCat event ID.
    """
    url = DETAIL_TEMPLATE.replace("[EVENTID]", eventid)
    with urlopen(url, timeout=60) as fh:
        data = fh.read().decode("utf8")
    jdict = json.loads(data)
    return jdict


def get_bytes(url):
    """Get simple bytes from a url."""
    with urlopen(url, timeout=60) as fh:
        data = fh.read()

    return data
