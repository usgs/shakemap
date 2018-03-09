# stdlib imports
import os.path
from collections import OrderedDict
from datetime import timedelta
from io import StringIO 

# third party imports
from shakelib.utils.containers import ShakeMapInputContainer
from configobj import ConfigObj
from libcomcat.search import get_event_by_id, search
from libcomcat.classes import DetailEvent
import pandas as pd
from lxml import etree
import time

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths

TIMEWINDOW = 60 # number of seconds to search for event matching origin time
DEGWINDOW = 0.1 # distance in decimal degrees to search for event matching coordinates
MAGWINDOW = 0.2 # +/- magnitude threshold to search for matching events

required_columns = ['station','lat','lon','network']
channel_groups = [['[a-z]{2}e','[a-z]{2}n','[a-z]{2}z'],
                  ['h1','h2','z'],
                  ['unk']]
pgm_cols = ['pga','pgv','psa03','psa10','psa30']
optional = ['location','distance','reference','intensity','source']

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
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        datafile = os.path.join(datadir, 'shake_data.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # try to find the event by our event id
        try:
            detail = get_event_by_id(self._eventid)
        except:
            # load the container, get the basic event info
            container = ShakeMapInputContainer.load(datafile)
            origin = container.getRuptureObject()._origin
            lat = origin.lat
            lon = origin.lon
            depth = origin.depth
            mag = origin.magnitude
            etime = origin.time
            tstart = etime-timedelta(seconds=TIMEWINDOW)
            tend = etime+timedelta(seconds=TIMEWINDOW)
            minlat = lat - DEGWINDOW
            minlon = lon - DEGWINDOW
            maxlat = lat + DEGWINDOW
            maxlon = lon + DEGWINDOW
            minmag  = mag - MAGWINDOW
            maxmag  = mag + MAGWINDOW
            summary = search(starttime=tstart,
                             endtime=tend,
                             minlatitude=minlat,
                             maxlatitude=maxlat,
                             minlonitude=minlon,
                             maxlonitude=maxlon,
                             minmagnitude=minmag,
                             maxmagnitude=maxmag)
            if not len(summary) or len(summary) > 1:
                self.logging.info('No single event found matching %s' % self._eventid)
                return

            detail = summary[0]

        
        dataframe,msg = _get_dyfi_dataframe(detail)
        if dataframe is None:
            self.logging.info(msg)
            return

        reference = 'USGS Did You Feel It? System'
        _dataframe_to_xml(dataframe,self._eventid,datadir,reference)

def _get_dyfi_dataframe(detail_or_url):
    if isinstance(detail_or_url,str):
        detail = DetailEvent(detail_or_url)
    else:
        detail = detail_or_url

    if not detail.hasProduct('dyfi'):
        msg = '%s has no DYFI product at this time.' % detail.url
        dataframe = None
        return (dataframe,msg)
    
    dyfi = detail.getProducts('dyfi')[0]
    text_data,url = dyfi.getContentBytes('cdi_zip.txt')
    text_data = text_data.decode('utf-8')
    lines = text_data.split('\n')
    columns = lines[0].split(':')[1].split(',')
    columns = [col.strip() for col in columns]
    fileio = StringIO(text_data)
    df = pd.read_csv(fileio,skiprows=1,names=columns)

    # we need to rename the dataframe columns to match what the function
    # below expects: station,lat,lon,network,intensity
    df = df.rename(index=str,columns={'ZIP/Location':'station',
                                      'City':'location',
                                      'CDI':'intensity',
                                      'Latitude':'lat',
                                      'Longitude':'lon',
                                      'Hypocentral distance':'distance'})

    df = df.drop(['No. of responses','Suspect?','State'],axis=1)
    df['network'] = 'DYFI'
    df['source'] = source="USGS (Did You Feel It?)"
    df = df.astype({'station':str})
    return (df,'')
        
def _dataframe_to_xml(df,eventid,dir,reference=None):
    """Write a dataframe to ShakeMap XML format.
    
    Args:
        df (DataFrame): Pandas dataframe, as described in read_excel.
        eventid (str): Event ID string.
        dir (str): Path to directory where XML file should be written.
    Returns:
        str: Path to output XML file.
    """
    #######################################################################
    # TODO: This function is a copy of one from amptools.  Resolve this
    # in some intelligent way.
    #######################################################################
    if hasattr(df.columns,'levels'):
        top_headers = df.columns.levels[0]
        channels = (set(top_headers) - set(required_columns)) - set(optional)
    else:
        channels = []
    root = etree.Element('shakemap-data',code_version="3.5",map_version="3")

    create_time = int(time.time())
    stationlist = etree.SubElement(root,'stationlist',created='%i' % create_time)
    if reference is not None:
        stationlist.attrib['reference'] = reference
    
    for idx,row in df.iterrows():
        station = etree.SubElement(stationlist,'station')

        tmprow = row.copy()
        if isinstance(tmprow.index,pd.core.indexes.multi.MultiIndex):
            tmprow.index = tmprow.index.droplevel(1)
        
        # assign required columns
        stationcode = tmprow['station'].strip()
        netid = tmprow['network'].strip()
        if not stationcode.startswith(netid):
            stationcode = '%s.%s' % (netid,stationcode)

        station.attrib['code'] = stationcode
        station.attrib['lat'] = '%.4f' % tmprow['lat']
        station.attrib['lon'] = '%.4f' % tmprow['lon']

        # assign optional columns
        if 'location' in tmprow:
            station.attrib['name'] = tmprow['location'].strip()
        if 'network' in tmprow:
            station.attrib['netid'] = tmprow['network'].strip()
        if 'distance' in tmprow:
            station.attrib['dist'] = '%.1f' % tmprow['distance']
        if 'intensity' in tmprow:
            station.attrib['intensity'] = '%.1f' % tmprow['intensity']
        if 'source' in tmprow:
            station.attrib['source'] = tmprow['source'].strip()

        # sort channels by N,E,Z or H1,H2,Z
        channels = sorted(list(channels))
            
        for channel in channels:
            component = etree.SubElement(station,'comp')
            component.attrib['name'] = channel.upper()

            # create sub elements out of any of the PGMs
            for pgm in ['pga','pgv','psa03','psa10','psa30']:
                if pgm not in row[channel] or np.isnan(row[channel][pgm]):
                    continue
                if pgm in row[channel]:
                    pgm_el = etree.SubElement(component,pgm)
                    pgm_el.attrib['flag'] = '0'
                    pgm_el.attrib['value'] = '%.4f' % row[channel][pgm]

    
    outfile = os.path.join(dir,'%s_dat.xml' % eventid)
    tree = etree.ElementTree(root)
    tree.write(outfile,pretty_print=True)

    return outfile
