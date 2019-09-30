# stdlib imports
import os.path
from io import StringIO
import json

# third party imports
from libcomcat.search import get_event_by_id
from libcomcat.classes import DetailEvent
from impactutils.io.table import dataframe_to_xml
import pandas as pd
import numpy as np

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths

# Get rid of stupid pandas warning
pd.options.mode.chained_assignment = None

# number of seconds to search for event matching origin time
TIMEWINDOW = 60
# distance in decimal degrees to search for event matching coordinates
DEGWINDOW = 0.1
# +/- magnitude threshold to search for matching events
MAGWINDOW = 0.2

required_columns = ['station', 'lat', 'lon', 'network']
channel_groups = [['[a-z]{2}e', '[a-z]{2}n', '[a-z]{2}z'],
                  ['h1', 'h2', 'z'],
                  ['unk']]
pgm_cols = ['pga', 'pgv', 'psa03', 'psa10', 'psa30']
optional = ['location', 'distance', 'reference', 'intensity', 'source']

# what are the DYFI columns and what do we rename them to?
DYFI_COLUMNS_REPLACE = {
    'Geocoded box': 'station',
    'CDI': 'intensity',
    'Latitude': 'lat',
    'Longitude': 'lon',
    'No. of responses': 'nresp',
    'Hypocentral distance': 'distance'
}

OLD_DYFI_COLUMNS_REPLACE = {
    'ZIP/Location': 'station',
    'CDI': 'intensity',
    'Latitude': 'lat',
    'Longitude': 'lon',
    'No. of responses': 'nresp',
    'Epicentral distance': 'distance'
}

MIN_RESPONSES = 3  # minimum number of DYFI responses per grid


class DYFIModule(CoreModule):
    """
    dyfi -- Search ComCat for DYFI data and turn it into a ShakeMap data file.
    """

    command_name = 'dyfi'

    def execute(self):
        """
        Write info.json metadata file.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        _, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        if not os.path.isdir(datadir):
            os.makedirs(datadir)

        # try to find the event by our event id
        try:
            detail = get_event_by_id(self._eventid)
            dataframe, msg = _get_dyfi_dataframe(detail)
        except Exception as e:
            fmt = 'Could not retrieve DYFI data for %s - error "%s"'
            self.logger.warning(fmt % (self._eventid, str(e)))
            return

        if dataframe is None:
            self.logger.info(msg)
            return

        reference = 'USGS Did You Feel It? System'
        xmlfile = os.path.join(datadir, 'dyfi_dat.xml')
        dataframe_to_xml(dataframe, xmlfile, reference)
        self.logger.info('Wrote %i DYFI records to %s' %
                         (len(dataframe), xmlfile))


def _get_dyfi_dataframe(detail_or_url, inputfile=None):

    if inputfile:
        with open(inputfile,'rb') as f:
            rawdata = f.read()
        if 'json' in inputfile:
            df = _parse_geocoded_json(rawdata)
        else:
            df = _parse_geocoded_csv(rawdata)
        if df is None:
            msg = 'Could not read file %s' % inputfile

    else:
        if isinstance(detail_or_url, str):
            detail = DetailEvent(detail_or_url)
        else:
            detail = detail_or_url

        df, msg = _parse_dyfi_detail(detail)

    if df is None:
        return None, msg

    df['netid'] = 'DYFI'
    df['source'] = "USGS (Did You Feel It?)"
    df.columns = df.columns.str.upper()

    return (df, '')


def _parse_dyfi_detail(detail):
    
    if not detail.hasProduct('dyfi'):
        msg = '%s has no DYFI product at this time.' % detail.url
        dataframe = None
        return (dataframe, msg)

    dyfi = detail.getProducts('dyfi')[0]

    # search the dyfi product, see which of the geocoded
    # files (1km or 10km) it has.  We're going to select the data from
    # whichever of the two has more entries with >= 3 responses,
    # preferring 1km if there is a tie.
    df_10k = pd.DataFrame({'a': []})
    df_1k = pd.DataFrame({'a': []})

    # get 10km data set, if exists
    if len(dyfi.getContentsMatching('dyfi_geo_10km.geojson')):
        bytes_10k, _ = dyfi.getContentBytes('dyfi_geo_10km.geojson')
        df_10k = _parse_geocoded_json(bytes_10k)

    # get 1km data set, if exists
    if len(dyfi.getContentsMatching('dyfi_geo_1km.geojson')):
        bytes_1k, _ = dyfi.getContentBytes('dyfi_geo_1km.geojson')
        df_1k = _parse_geocoded_json(bytes_1k)

    if len(df_1k) >= len(df_10k):
        df = df_1k
    else:
        df = df_10k

    if not len(df):
        # try to get a text file data set
        if not len(dyfi.getContentsMatching('cdi_geo.txt')):
            return (None, 'No geocoded datasets are available for this event.')

        bytes_geo, _ = dyfi.getContentBytes('cdi_geo.txt')
        df = _parse_geocoded_csv(bytes_geo)
        
    return df, ''


def _parse_geocoded_csv(bytes_data):
    # the dataframe we want has columns:
    # 'intensity', 'distance', 'lat', 'lon', 'station', 'nresp'
    # the cdi geo file has:
    # Geocoded box, CDI, No. of responses, Hypocentral distance,
    # Latitude, Longitude, Suspect?, City, State

    # download the text file, turn it into a dataframe
    
    text_geo = bytes_data.decode('utf-8')
    lines = text_geo.split('\n')
    columns = lines[0].split(':')[1].split(',')
    columns = [col.strip() for col in columns]

    fileio = StringIO(text_geo)
    df = pd.read_csv(fileio, skiprows=1, names=columns)
    if 'ZIP/Location' in columns:
        df = df.rename(index=str, columns=OLD_DYFI_COLUMNS_REPLACE)
    else:
        df = df.rename(index=str, columns=DYFI_COLUMNS_REPLACE)
    df = df.drop(['Suspect?', 'City', 'State'], axis=1)
    df = df[df['nresp'] >= MIN_RESPONSES]

    return df


def _parse_geocoded_json(bytes_data):
    
    text_data = bytes_data.decode('utf-8')
    jdict = json.loads(text_data)
    if len(jdict['features']) == 0:
        return None
    prop_columns = list(jdict['features'][0]['properties'].keys())
    columns = ['lat', 'lon'] + prop_columns
    arrays = [[] for col in columns]
    df_dict = dict(zip(columns, arrays))
    for feature in jdict['features']:
        for column in prop_columns:
            if column == 'name':
                prop = feature['properties'][column]
                prop = prop[0:prop.find('<br>')]
            else:
                prop = feature['properties'][column]

            df_dict[column].append(prop)
        # the geojson defines a box, so let's grab the center point
        lons = [c[0] for c in feature['geometry']['coordinates'][0]]
        lats = [c[1] for c in feature['geometry']['coordinates'][0]]
        clon = np.mean(lons)
        clat = np.mean(lats)
        df_dict['lat'].append(clat)
        df_dict['lon'].append(clon)

    df = pd.DataFrame(df_dict)
    df = df.rename(index=str, columns={
        'cdi': 'intensity',
        'dist': 'distance',
        'name': 'station'
    })
    if df is not None:
        df = df[df['nresp'] >= MIN_RESPONSES]

    return df
