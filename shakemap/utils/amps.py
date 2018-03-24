import argparse
import sys
import sqlite3
import glob
import os.path
from xml.dom import minidom
from datetime import datetime
from collections import OrderedDict

# third party libraries
import numpy as np
import pandas as pd
from openquake.hazardlib.geo.geodetic import geodetic_distance

EVENT = OrderedDict([('id','INTEGER PRIMARY KEY'),
                     ('eventid','TEXT'),
                     ('directory','TEXT'),
                     ('directory_time','INTEGER'),
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

DBFILE = 'amps.db'

IMTS = ['acc','vel','sa','pga','pgv','pgd']
IMTDICT = {'acc':'pga',
           'vel':'pgv'}

TMIN = 30
TMAX = 180
DISTANCE = 500

class AmplitudeHandler(object):
    def __init__(self,install_path,data_path):
        self._data_path = data_path
        dbfile = os.path.join(install_path,'data',DBFILE)
        db_exists = os.path.isfile(dbfile)
        self._conn = sqlite3.connect(dbfile)
        self._cursor = self._conn.cursor()
        if not db_exists:
            for table,tdict in TABLES.items():
                createcmd = 'CREATE TABLE %s (' % table
                nuggets = []
                for column,ctype in tdict.items():
                    nuggets.append('%s %s' % (column,ctype))
                createcmd += ','.join(nuggets) + ')'
                self._cursor.execute(createcmd)

    def lock(self):
        self._cursor.execute('PRAGMA locking_mode = EXCLUSIVE')
        self._cursor.execute('BEGIN EXCLUSIVE')

    def unlock(self):
        self._conn.commit()

    def update(self):
        # update the event table with anything new on the file system
        equery = 'SELECT max(directory_time) FROM event'
        self._cursor.execute(equery)
        tlast = self._cursor.fetchone()[0] # unix timestamp with most recent known event
        if tlast is None: # event table is not populated
            tlast = 0
        event_dirs = os.listdir(self._data_path)
        nadded = 0
        for tdir in event_dirs:
            event_dir = os.path.join(self._data_path,tdir)
            if not os.path.isdir(event_dir):
                continue
            if not os.path.isdir(os.path.join(event_dir,'current')):
                continue
            dtime = os.path.getmtime(event_dir)
            if dtime > tlast:
                eventxml = os.path.join(event_dir,'current','event.xml')
                event = _read_event(eventxml)
                cols = '(eventid,directory,directory_time,time,lat,lon,depth,magnitude)'
                fmt = 'INSERT INTO event %s VALUES ("%s","%s",%i,%i,%.4f,%.4f,%.1f,%.1f)'
                einsert = fmt % (cols,
                                 event['id'],
                                 event_dir,
                                 dtime,
                                 event['time'],
                                 event['lat'],
                                 event['lon'],
                                 event['depth'],
                                 event['mag'])
                self._cursor.execute(einsert)
                nadded += 1

        # also delete any entries for events that are no longer on the file system
        equery2 = 'SELECT id,directory FROM event'
        self._cursor.execute(equery2)
        event_rows = self._cursor.fetchall()
        ndeleted = 0
        for event_row in event_rows:
            eid,efolder = event_row
            if not os.path.isdir(efolder):
                dquery = 'DELETE FROM event WHERE id=%i' % eid
                self._cursor.execute(dquery)
                ndeleted += 1
                
        return (nadded,ndeleted)

    def associateAll(self):
        # loop through events in database
        # find associated peak ground motions
        equery = 'SELECT directory,time,lat,lon FROM event'
        self._cursor.execute(equery)
        events = self._cursor.fetchall()
        nassociated = 0
        for event in events:
            efolder = event[0]
            eqtime = event[1]
            lat = event[2]
            lon = event[3]
            dataframe = self.associate(eqtime,lat,lon)
            if not len(dataframe):
                continue
            amptime = datetime.utcnow().strftime('%Y%m%d%H%M%S')
            outfile = os.path.join(efolder,'%s_amps.xlsx' % amptime)
            dataframe.to_excel(outfile,index=False)
            nassociated += 1
        
    def associate(self,eqtime,eqlat,eqlon):
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
        
        return dataframe

    def __del__(self):
        self.unlock()
        self._conn.close()
    
    def insert(self,xmlfile):
        fpath,fname = os.path.split(xmlfile)
        try:
            xmldata = open(xmlfile,'rb').read().decode('utf-8')
            root = minidom.parseString(xmldata)
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

def _read_event(xmlfile):
    edict = {}
    root = minidom.parse(xmlfile)
    eq = root.getElementsByTagName('earthquake')[0]
    year = int(eq.getAttribute('year'))
    month = int(eq.getAttribute('month'))
    day = int(eq.getAttribute('day'))
    hour = int(eq.getAttribute('hour'))
    minute = int(eq.getAttribute('minute'))
    second = int(eq.getAttribute('second'))
    etime = int(datetime(year,month,day,hour,minute,second).timestamp())
    lat = float(eq.getAttribute('lat'))
    lon = float(eq.getAttribute('lon'))
    depth = float(eq.getAttribute('depth'))
    mag = float(eq.getAttribute('mag'))
    eid = eq.getAttribute('id')
    net = eq.getAttribute('network')
    if not eid.startswith(net):
        eid = net+eid
    edict = {'id':eid,
             'time':etime,
             'lat':lat,
             'lon':lon,
             'depth':depth,
             'mag':mag}
    root.unlink()
    return edict
    
