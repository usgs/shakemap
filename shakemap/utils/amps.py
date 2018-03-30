import argparse
import sys
import sqlite3
import glob
import os.path
import re
from xml.dom import minidom
from datetime import datetime
from collections import OrderedDict

# third party libraries
import numpy as np
import pandas as pd
from openquake.hazardlib.geo.geodetic import geodetic_distance
from amptools.table import dataframe_to_xml

# define all of the tables as dictionaries
EVENT = OrderedDict([('id','INTEGER PRIMARY KEY'),
                     ('eventid','TEXT'),
                     ('time','INTEGER'),
                     ('lat','REAL'),
                     ('lon','REAL'),
                     ('depth','REAL'),
                     ('magnitude','REAL')])

STATION = OrderedDict([('id','INTEGER PRIMARY KEY'),
                       ('timestamp','INTEGER'),
                       ('lat','REAL'),
                       ('lon','REAL'),
                       ('network','TEXT'),
                       ('name','TEXT'),
                       ('code','TEXT')])

CHANNEL = OrderedDict([('id','INTEGER PRIMARY KEY'),
                       ('station_id','INTEGER'),
                       ('channel','TEXT')])

PGM = OrderedDict([('id','INTEGER PRIMARY KEY'),
                   ('channel_id','INTEGER'),
                   ('imt','TEXT'),
                   ('value','REAL')])

TABLES = {'event':EVENT,
          'station':STATION,
          'channel':CHANNEL,
          'pgm':PGM}

# 2018-03-21T21:19:16.625Z
TIMEFMT = '%Y-%m-%dT%H:%M:%S'

# database file name
DBFILE = 'amps.db'

IMTS = ['acc','vel','sa','pga','pgv']
# sometimes (sigh) pga/pgv labeled as acc/vel
IMTDICT = {'acc':'pga',
           'vel':'pgv'}

# association algorithm - any peak with:
# time > origin - TMIN and time < origin + TMAX
# AND
# distance < DISTANCE
TMIN = 30
TMAX = 180
DISTANCE = 500

class AmplitudeHandler(object):
    """Store and associate strong motion peak amplitudes with earthquake events.

    """
    MAX_CONNECT_ATTEMPTS = 3
    CONNECT_TIMEOUT = 5 # number of seconds to try getting past lock
    def __init__(self,install_path,data_path):
        """Instantiate amplitude handler with ShakeMap profile paths.

        """
        self._data_path = data_path
        self._dbfile = os.path.join(install_path,'data',DBFILE)
        db_exists = os.path.isfile(self._dbfile)
        self._connection = None
        self._cursor = None
        self._connect()
        if not db_exists:
            for table,tdict in TABLES.items():
                createcmd = 'CREATE TABLE %s (' % table
                nuggets = []
                for column,ctype in tdict.items():
                    nuggets.append('%s %s' % (column,ctype))
                createcmd += ','.join(nuggets) + ')'
                self._cursor.execute(createcmd)
        self._disconnect()

    def _connect(self):
        for i in range(0,3):
            try:
                self._connection = sqlite3.connect(self._dbfile)
                self._cursor = self._connection.cursor()
                break
            except:
                continue
        if self._connection is None:
            raise Exception('Could not connect to %s' % self._dbfile)

        self._lock()

    def _disconnect(self):
        self._unlock()
        self._cursor.close()
        self._connection.close()
        self._connection = None
        self._cursor = None

    def _lock(self):
        """Lock the sqlite database to prevent other processes from opening the file.

        """
        self._cursor.execute('PRAGMA locking_mode = EXCLUSIVE')
        self._cursor.execute('BEGIN EXCLUSIVE')

    def _unlock(self):
        """UnLock the sqlite database and commit the changes made since locking.

        """
        self._connection.commit()


    def insertEvent(self,event):
        """Insert an event into the database.
        
        An directory with name of event['id'] should exist in data_path.

        Args:
            event (dict): Dictionary containing fields:
                          - id: Event ID (i.e., us2008abcd).
                          - time: Origin time in UTC (datetime).
                          - lat: Origin latitude (dd).
                          - lon: Origin longitude (dd).
                          - depth: Origin depth (km).
                          - mag: Earthquake magnitude.
        """
        self._connect()

        cols = '(eventid,time,lat,lon,depth,magnitude)'
        fmt = 'INSERT INTO event %s VALUES ("%s",%i,%.4f,%.4f,%.1f,%.1f)'
        einsert = fmt % (cols,
                         event['id'],
                         event['time'].timestamp(),
                         event['lat'],
                         event['lon'],
                         event['depth'],
                         event['mag'])
        self._cursor.execute(einsert)
        
        self._disconnect()

    def associateAll(self,write_data=True):
        """Associate peak ground motions with appropriate events, write station XML to file system.
        
        Ground motion records associated with events will be deleted from the database.

        Args:
            write_data (bool): Control whether associated data is written to directory.
                               Used for testing.
        
        Returns:
            int: Number of events with associated ground motions.
        """
        self._connect()
        equery = 'SELECT eventid,time,lat,lon FROM event'
        self._cursor.execute(equery)
        events = self._cursor.fetchall()
        self._disconnect()
        
        nassociated = 0
        for event in events:
            eventid = event[0]
            eqtime = event[1]
            lat = event[2]
            lon = event[3]
            dataframe = self.associate(eqtime,lat,lon)
            if not len(dataframe):
                continue
            amptime = datetime.utcnow().strftime('%Y%m%d%H%M%S')
            data_folder = os.path.join(self._data_path,eventid,'current')
            if write_data:
                if not os.path.isdir(data_folder):
                    os.makedirs(data_folder)
                xmlfile = os.path.join(data_folder,'unassoc_%s_dat.xml' % amptime)
                dataframe_to_xml(dataframe,xmlfile)
                
            nassociated += 1

        return nassociated
        
    def associate(self,eqtime,eqlat,eqlon):
        """Find peak ground motion records associated with input event info.

        Ground motion records associated with input event are deleted from the
        database.

        Args:
            eqtime (int): Unix timestamp of earthquake origin.
            eqlat (float): Latitude of earthquake origin.
            eqlon (float): Longitude of earthquake origin.
        Returns:
            DataFrame: Pandas dataframe containing peak ground motions.  Columns:
                       - station: Station code
                       - channel: Channel (HHE,HHN, etc.)
                       - imt: Intensity measure type (pga,pgv, etc.)
                       - value: IMT value.
                       - lat: Station latitude.
                       - lon: Station longitude.
                       - network: Station contributing network.
                       - location: String describing station location.
                       - distance: Distance (km) from station to origin.
        """
        self._connect()
        columns = ['station', 'channel','imt','value',
               'lat', 'lon', 'network','location', 'distance']
        dataframe = pd.DataFrame(columns=columns)
        conditions = ['timestamp > %i' % (eqtime-TMIN),
                      'timestamp < %i' % (eqtime+TMAX)]
        cond_str = ' AND '.join(conditions)
        time_cmd = 'SELECT id,timestamp,lat,lon FROM station WHERE %s ORDER by timestamp DESC' % cond_str
        self._cursor.execute(time_cmd)
        # numpy array of id, timestamp, lat, lon
        eqdata = np.array(self._cursor.fetchall())

        if not len(eqdata):
            return dataframe
        
        dist = geodetic_distance(eqlon,eqlat,eqdata[:,3],eqdata[:,2])

        inear = np.where(dist < DISTANCE)[0]
        eqdata = eqdata[inear]
        dist = dist[inear]
        nrows,ncols = eqdata.shape

        stations = []

        # make arrays of rows to delete after we've done association
        station_ids = []
        channel_ids = []
        pgm_ids = []
        
        for idx in range(0,nrows):
            station_id = int(eqdata[idx,0])
            station_ids.append(station_id)
            fmt = 'SELECT network,name,code FROM station WHERE id=%i'
            squery = fmt % station_id
            self._cursor.execute(squery)
            station = self._cursor.fetchone()
            code = station[2]
            if code in stations:
                continue
            stations.append(code)
            lat = eqdata[idx,2]
            lon = eqdata[idx,3]
            distance = dist[idx]
            network = station[0]
            name = station[1]

            # find all the channels associated with this station
            fmt = 'SELECT id,channel FROM channel WHERE station_id=%i'
            cquery = fmt % station_id
            self._cursor.execute(cquery)
            channels = self._cursor.fetchall()
            for channel in channels:
                channel_id = channel[0]
                channel_ids.append(channel_id)
                channel_name = channel[1]

                # find all of the pgms associated with that channel
                fmt = 'SELECT id,imt,value FROM pgm WHERE channel_id=%i'
                pquery = fmt % channel_id
                self._cursor.execute(pquery)
                pgms = self._cursor.fetchall()
                for pgm in pgms:
                    pgm_id = pgm[0]
                    pgm_ids.append(pgm_id)
                    imt = pgm[1]
                    value = pgm[2]
                    data_row = pd.Series(index=columns)
                    data_row['station'] = code
                    data_row['channel'] = channel_name
                    data_row['imt'] = imt
                    data_row['value'] = value
                    data_row['lat'] = lat
                    data_row['lon'] = lon
                    data_row['network'] = network
                    data_row['location'] = name
                    data_row['distance'] = distance
                    dataframe = dataframe.append(data_row,ignore_index=True)

        # clean up rows that have been associated
        for pgm_id in pgm_ids:
            delete_cmd = 'DELETE FROM pgm WHERE id=%i' % pgm_id
            self._cursor.execute(delete_cmd)
            
        for channel_id in channel_ids:
            delete_cmd = 'DELETE FROM channel WHERE id=%i' % channel_id
            self._cursor.execute(delete_cmd)

        for station_id in station_ids:
            delete_cmd = 'DELETE FROM station WHERE id=%i' % station_id
            self._cursor.execute(delete_cmd)

        # # ensure that the station field is a string
        dataframe['station']= dataframe['station'].astype(int)
        dataframe['station']= dataframe['station'].astype(str)

        self._disconnect()
        
        return dataframe

    def __del__(self):
        """Destructor.

        """
        if self._connection is not None:
            self._disconnect()
    
    def insertAmps(self,xmlfile):
        """Insert data from amps file into database.
        
        Args:
            xmlfile (str): XML file containing peak ground motion data.
            
        """
        self._connect()
        fpath,fname = os.path.split(xmlfile)
        try:
            xmlstr = open(xmlfile,'r').read()
            # sometimes these records have non-ascii bytes in them
            newxmlstr = re.sub(r'[^\x00-\x7F]+',' ', xmlstr)
            # newxmlstr = _invalid_xml_remove(xmlstr)
            newxmlstr = newxmlstr.encode('utf-8',errors='xmlcharrefreplace')
            root = minidom.parseString(newxmlstr)
        except Exception as e:
            raise Exception('Could not parse %s, due to error "%s"' % (xmlfile,str(e)))

        try:
            amps = root.getElementsByTagName('amplitudes')[0]
        except:
            raise Exception('%s does not appear to be an amplitude XML file.' % xmlfile)
        agency = amps.getAttribute('agency')
        record = amps.getElementsByTagName('record')[0]
        timing = record.getElementsByTagName('timing')[0]
        reference = timing.getElementsByTagName('reference')[0]
        has_pgm = False
        time_dict = {}
        for child in reference.childNodes:
            if child.nodeName == 'PGMTime':
                has_pgm = True
            if child.nodeName == 'year':
                time_dict['year'] = int(child.getAttribute('value'))
            if child.nodeName == 'month':
                time_dict['month'] = int(child.getAttribute('value'))
            if child.nodeName == 'day':
                time_dict['day'] = int(child.getAttribute('value'))
            if child.nodeName == 'hour':
                time_dict['hour'] = int(child.getAttribute('value'))
            if child.nodeName == 'minute':
                time_dict['minute'] = int(child.getAttribute('value'))
            if child.nodeName == 'second':
                time_dict['second'] = int(child.getAttribute('value'))
            if child.nodeName == 'msec':
                time_dict['msec'] = int(child.getAttribute('value'))
        if has_pgm:
            pgmtime_str = reference.getElementsByTagName('PGMTime')[0].firstChild.nodeValue
            pgmtime = int(datetime.strptime(pgmtime_str[0:19],TIMEFMT).timestamp())
        else:
            if not len(time_dict):
                print('No time data for file %s' % fname)
                return
            pgmtime = int(datetime(time_dict['year'],
                               time_dict['month'],
                               time_dict['day'],
                               time_dict['hour'],
                               time_dict['minute'],
                               time_dict['second']).timestamp())

        # there are often multiple stations per file, but they're all duplicates
        # of each other, so just grab the information from the first one
        station = record.getElementsByTagName('station')[0]
        lat = float(station.getAttribute('lat'))
        lon = float(station.getAttribute('lon'))
        code = station.getAttribute('code')
        name = station.getAttribute('name')
        if station.hasAttribute('net'):
            network = station.getAttribute('net')
        else:
            network = agency

        cols = '(timestamp,lat,lon,name,code,network)'
        fmt = 'INSERT INTO station %s VALUES (%i,%.4f,%.4f,"%s","%s","%s")'
        insert_cmd = fmt % (cols,pgmtime,lat,lon,name,code,network)
        self._cursor.execute(insert_cmd)
        sid = self._cursor.lastrowid

        # loop over components
        for channel in record.getElementsByTagName('component'):
            cname = channel.getAttribute('name')
            cols = '(station_id,channel)'
            fmt = 'INSERT INTO channel %s VALUES (%i,"%s")'
            insert_cmd = fmt % (cols,sid,cname)
            self._cursor.execute(insert_cmd)
            cid = self._cursor.lastrowid

            # loop over imts in channel
            for pgm in channel.childNodes:
                imt = pgm.nodeName
                if imt not in IMTS:
                    continue
                if imt == 'sa':
                    imt = 'p'+imt+pgm.getAttribute('period').replace('.','')
                if imt in IMTDICT:
                    imt = IMTDICT[imt]
                value = float(pgm.getAttribute('value'))
                cols = '(channel_id,imt,value)'
                fmt = 'INSERT INTO pgm %s VALUES (%i,"%s",%.4f)'
                insert_cmd = fmt % (cols,cid,imt,value)
                self._cursor.execute(insert_cmd)
                pid = self._cursor.lastrowid

        root.unlink()
        self._disconnect()

    def cleanAmps(self,threshold=30):
        """Clean out amplitude data that is older than the threshold number of days.

        Args:
            threshold (int): Maximum age in days of amplitude data in the database.
        Returns:
            int: Number of stations deleted.
        """
        tnow = int(datetime.utcnow().timestamp())
        self._connect()
        squery = 'SELECT id FROM station WHERE timestamp < %i' % tnow
        self._cursor.execute(squery)
        station_ids = self._cursor.fetchall()
        for station_id in station_ids:
            cquery = 'SELECT id FROM channel WHERE station_id=%i' % station_id
            self._cursor.execute(squery)
            channel_ids = self._cursor.fetchall()
            for channel_id in channel_ids:
                pquery = 'SELECT id FROM pgm WHERE channel_id=%i' % channel_id
                self._cursor.execute(pquery)
                pmg_ids = self._cursor.fetchall()
                for pgm_id in pgm_ids:
                    pdelete = 'DELETE FROM pmg WHERE id=%i' % pgm_id
                    self.cursor.execute(pdelete)
                cdelete = 'DELETE FROM channel WHERE id=%i' % channel_id
            sdelete = 'DELETE FROM station WHERE id=%i' % station_id
            self.cursor.execute(sdelete)

        self._disconnect()
        return len(station_ids)

    def cleanEvents(self,threshold=365):
        """Clean out event data that is older than the threshold number of days.

        Args:
            threshold (int): Maximum age in days of events in the database.
        Returns:
            int: Number of events deleted.
        """
        self._connect()
        tnow = int(datetime.utcnow().timestamp())
        equery = 'DELETE FROM event WHERE time < %i' % tnow
        self._cursor.execute(equery)
        nevents = self._cursor.rowcount
        self._disconnect()
        return nevents

    def getStats(self):
        """Get summary statistics about the database.

        Returns:
            dict: Fields: 
                  - events Number of events in database.
                  - stations Number of stations in database.
                  - channels Number of unique channels in database.
                  - pgms Number of unique pgms in database.
                  - event_min: Datetime of earliest event in database.
                  - event_max: Datetime of most recent event in database.
                  - station_min: Datetime of earliest amplitude data in database.
                  - station_max: Datetime of most recent amplitude data in database.

        """
        self._connect()
        results = {}

        # event stuff
        equery = 'SELECT count(*),min(time),max(time) FROM event'
        self._cursor.execute(equery)
        row = self._cursor.fetchone()
        results['events'] = row[0]
        if row[0] == 0:
            results['event_min'] = None
            results['event_max'] = None
        else:
            results['event_min'] = datetime.utcfromtimestamp(row[1])
            results['event_max'] = datetime.utcfromtimestamp(row[2])

        # station stuff
        squery = 'SELECT count(*), min(timestamp), max(timestamp) FROM station'
        self._cursor.execute(squery)
        row = self._cursor.fetchone()
        results['stations'] = row[0]
        if row[0] == 0:
            results['station_min'] = None
            results['station_max'] = None
        else:
            results['station_min'] = datetime.utcfromtimestamp(row[1])
            results['station_max'] = datetime.utcfromtimestamp(row[2])

        # channels
        cquery = 'SELECT count(*) FROM channel'
        self._cursor.execute(cquery)
        row = self._cursor.fetchone()
        results['channels'] = row[0]

        # pgms
        pquery = 'SELECT count(*) FROM pgm'
        self._cursor.execute(pquery)
        row = self._cursor.fetchone()
        results['pgms'] = row[0]

        self._disconnect()
        return results

def _invalid_xml_remove(c):
    #http://stackoverflow.com/questions/1707890/fast-way-to-filter-illegal-xml-unicode-chars-in-python
    illegal_unichrs = [ (0x00, 0x08), (0x0B, 0x1F), (0x7F, 0x84), (0x86, 0x9F),
                    (0xD800, 0xDFFF), (0xFDD0, 0xFDDF), (0xFFFE, 0xFFFF),
                    (0x1FFFE, 0x1FFFF), (0x2FFFE, 0x2FFFF), (0x3FFFE, 0x3FFFF),
                    (0x4FFFE, 0x4FFFF), (0x5FFFE, 0x5FFFF), (0x6FFFE, 0x6FFFF),
                    (0x7FFFE, 0x7FFFF), (0x8FFFE, 0x8FFFF), (0x9FFFE, 0x9FFFF),
                    (0xAFFFE, 0xAFFFF), (0xBFFFE, 0xBFFFF), (0xCFFFE, 0xCFFFF),
                    (0xDFFFE, 0xDFFFF), (0xEFFFE, 0xEFFFF), (0xFFFFE, 0xFFFFF),
                    (0x10FFFE, 0x10FFFF) ]

    illegal_ranges = ["%s-%s" % (unichr(low), unichr(high)) 
                  for (low, high) in illegal_unichrs 
                  if low < sys.maxunicode]

    illegal_xml_re = re.compile(u'[%s]' % u''.join(illegal_ranges))
    if illegal_xml_re.search(c) is not None:
        #Replace with space
        return ' '
    else:
        return c
