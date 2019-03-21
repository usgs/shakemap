# import sys
import sqlite3
import os.path
import re
from xml.dom import minidom
from datetime import datetime, timezone, timedelta
from collections import OrderedDict
import time
import defusedxml.cElementTree as dET
import xml.etree.cElementTree as ET
import json
from itertools import zip_longest

# third party libraries
import numpy as np
from openquake.hazardlib.geo.geodetic import geodetic_distance

# local libraries
from shakelib.rupture import constants

# define all of the tables as dictionaries
EVENT = OrderedDict([('id', 'INTEGER PRIMARY KEY'),
                     ('eventid', 'TEXT UNIQUE'),
                     ('netid', 'TEXT'),
                     ('network', 'TEXT'),
                     ('time', 'INTEGER'),
                     ('lat', 'REAL'),
                     ('lon', 'REAL'),
                     ('depth', 'REAL'),
                     ('magnitude', 'REAL'),
                     ('locstring', 'TEXT'),
                     ('repeats', 'TEXT'),
                     ('lastrun', 'INTEGER')])

STATION = OrderedDict([('id', 'INTEGER PRIMARY KEY'),
                       ('timestamp', 'INTEGER'),
                       ('lat', 'REAL'),
                       ('lon', 'REAL'),
                       ('network', 'TEXT'),
                       ('name', 'TEXT'),
                       ('code', 'TEXT')])

CHANNEL = OrderedDict([('id', 'INTEGER PRIMARY KEY'),
                       ('station_id',
                        'INTEGER REFERENCES station(id) ON DELETE CASCADE'),
                       ('channel', 'TEXT'),
                       ('loc', 'TEXT')])

PGM = OrderedDict([('id', 'INTEGER PRIMARY KEY'),
                   ('channel_id',
                    'INTEGER REFERENCES channel(id) ON DELETE CASCADE'),
                   ('imt', 'TEXT'),
                   ('value', 'REAL')])

TABLES = {'event': EVENT,
          'station': STATION,
          'channel': CHANNEL,
          'pgm': PGM}

# database file name
DBFILE = 'amps.db'

IMTS = ['acc', 'vel', 'sa', 'pga', 'pgv']
# sometimes (sigh) pga/pgv labeled as acc/vel
IMTDICT = {'acc': 'pga',
           'vel': 'pgv'}

# association algorithm - any peak with:
# time > origin - TMIN and time < origin + TMAX
# AND
# distance < DISTANCE
TMIN = 60
TMAX = 180
DISTANCE = 500

# SQLite has a limit (999) on the number of variables in
# a query; we set our threshold somewhat lower than that for
# safety.
MAX_VARS = 200


class AmplitudeHandler(object):
    """Store and associate strong motion peak amplitudes with
    earthquake events.
    """

    def __init__(self, install_path, data_path):
        """Instantiate amplitude handler with ShakeMap profile paths.

        """
        self._data_path = data_path
        self._dbfile = os.path.join(install_path, 'data', DBFILE)
        db_exists = os.path.isfile(self._dbfile)
        self._connect()
        if not db_exists:
            for table, tdict in TABLES.items():
                createcmd = 'CREATE TABLE %s (' % table
                nuggets = []
                for column, ctype in tdict.items():
                    nuggets.append('%s %s' % (column, ctype))
                createcmd += ','.join(nuggets) + ')'
                self._cursor.execute(createcmd)
            self._cursor.execute('CREATE INDEX station_index ON '
                                 'channel(station_id)')
            self._cursor.execute('CREATE INDEX channel_index ON '
                                 'pgm(channel_id)')
            self._cursor.execute('CREATE INDEX eventid_index ON '
                                 'event(eventid)')
            self._cursor.execute('CREATE INDEX stacode_index ON '
                                 'station(code)')
            self._cursor.execute('CREATE INDEX stanet_index ON '
                                 'station(network)')
            self._cursor.execute('PRAGMA journal_mode = WAL')

    def _connect(self):
        self._connection = sqlite3.connect(self._dbfile, timeout=15)
        if self._connection is None:
            raise RuntimeError('Could not connect to %s' % self._dbfile)
        self._connection.isolation_level = 'EXCLUSIVE'
        self._cursor = self._connection.cursor()
        self._cursor.execute('PRAGMA foreign_keys = ON')
        self._cursor.execute('PRAGMA journal_mode = WAL')

    def _disconnect(self):
        self.commit()
        self._cursor.close()
        self._connection.close()
        self._connection = None
        self._cursor = None

    def commit(self):
        """Commit any operations to the database.
        """
        self._connection.commit()

    def insertEvent(self, event, update=False):
        """Insert an event into the database.

        A directory with name of event['id'] should exist in data_path.

        Args:
            event (dict): Dictionary containing fields:
                          - id: Event ID (i.e., us2008abcd).
                          - netid: Network code (i.e., us).
                          - network: Network name (i.e., "USGS Network").
                          - time: Origin time in UTC (datetime).
                          - lat: Origin latitude (dd).
                          - lon: Origin longitude (dd).
                          - depth: Origin depth (km).
                          - mag: Earthquake magnitude.
                          - locstring: Location string (i.e. '2 mi SE of Reno')
                          - repeats: A list of repeat times (optional)
                          - lastrun: Timestamp of the last run of the event.
                                     (optional)
            update (bool): Update an existing event with new info (True) or
                           insert a new event (False)

        Returns:
            nothing: Nothing.
        """
        cols = [x for x in EVENT.keys() if x != 'id']
        if update:
            # This makes a string like 'eventid = ?, netid = ?, ...'
            einsert = ('UPDATE event SET '
                       + ', '.join([' = '.join(x) for x in zip_longest(
                           cols, [], fillvalue='?')])
                       + ' WHERE eventid = "'
                       + str(event['id'])
                       + '"')
        else:
            einsert = ('INSERT INTO event ('
                       + ', '.join(cols)
                       + ') VALUES ('
                       + ', '.join('?' * len(cols))
                       + ')')
        if 'network' in event:
            network = event['network']
        else:
            network = ''
        if 'repeats' in event and event['repeats'] and \
           len(event['repeats']) > 0:
            repeats = json.dumps(event['repeats'])
        else:
            repeats = None
        if 'lastrun' in event:
            lastrun = event['lastrun']
        else:
            lastrun = int(time.time())
        self._cursor.execute(einsert, (event['id'],
                                       event['netid'],
                                       network,
                                       timestr_to_timestamp(event['time']),
                                       event['lat'],
                                       event['lon'],
                                       event['depth'],
                                       event['mag'],
                                       event['locstring'],
                                       repeats,
                                       lastrun))

        self.commit()
        return

    def getEvent(self, eventid):
        """Return the event parameters for the specified event.

        Args:
            eventid (str): The id of the event to query

        Returns:
            dictionary: A dictionary of the columns of the table and
            their values for the the event; a None is
            returned if the event is not in the database.
        """
        query = 'SELECT * FROM event WHERE eventid = ?'
        self._cursor.execute(query, (eventid,))
        row = self._cursor.fetchone()
        if row is None:
            return None
        cols = [col[0] for col in self._cursor.description]
        event = dict(zip(cols, row))
        #
        # Deal with differences between the database column names
        # and the event keys
        #
        event['id'] = event['eventid']
        del event['eventid']
        event['mag'] = event['magnitude']
        del event['magnitude']
        event['time'] = datetime.fromtimestamp(event['time'], timezone.utc).\
            strftime(constants.TIMEFMT)
        if event['repeats']:
            event['repeats'] = json.loads(event['repeats'])
        return event

    def deleteEvent(self, eventid):
        """Delete the event from the database.

        Args:
            eventid (str): The id of the event to delete

        Returns:
            nothing: Nothing.
        """
        query = 'DELETE FROM event WHERE eventid = ?'
        self._cursor.execute(query, (eventid,))
        self.commit()
        return

    def getRepeats(self):
        """Return all the rows from the event table where the 'repeats' column
        is not NULL.

        Args:
            none

        Returns:
            (list): List of tuples of (eventid, origin_time, [repeats]).
        """
        query = "SELECT eventid, time, repeats FROM event WHERE "\
                "repeats IS NOT NULL"
        self._cursor.execute(query)
        repeats = self._cursor.fetchall()
        replist = []
        for repeat in repeats:
            rep = list(repeat)
            rep[2] = json.loads(rep[2])
            replist.append(rep)
        return replist

    def associateAll(self, pretty_print=False):
        """Associate peak ground motions with appropriate events,
        write station XML to file system.

        Ground motion records associated with events will be deleted
        from the database.

        Args:
            pretty_print (bool): Writes more human-readable XML, but is
                                 slower and writes larger files. False
                                 by default.

        Returns:
            list: The event IDs of the events for which associated data
            were found.
        """
        equery = 'SELECT eventid, time, lat, lon FROM event'
        self._cursor.execute(equery)
        events = self._cursor.fetchall()

        associated = []
        for event in events:
            eventid = event[0]
            eqtime = event[1]
            lat = event[2]
            lon = event[3]

            data_list = self.associate(eqtime, lat, lon)

            if len(data_list) == 0:
                continue

            self.writeXML(data_list, eventid, pretty_print)

            associated.append(eventid)

        return associated

    def associateOne(self, eventid, pretty_print=False):
        """Associate peak ground motions with the specified event,
        write station XML to file system.

        Ground motion records associated with events will be deleted
        from the database.

        Args:
            eventid (str): The event ID of the event to associate
            pretty_print (bool): Writes more human-readable XML, but is
                                 slower and writes larger files. False
                                 by default.

        Returns:
            int: The number of amps associated with the specified event.
            -1 is returned if the event is not found in the database.
        """
        equery = 'SELECT time, lat, lon FROM event where eventid = ?'
        self._cursor.execute(equery, (eventid,))
        event = self._cursor.fetchone()

        if event is None:
            return -1

        data_list = self.associate(event[0], event[1], event[2])

        namps = len(data_list)

        if namps == 0:
            return 0

        self.writeXML(data_list, eventid, pretty_print)

        return namps

    def associate(self, eqtime, eqlat, eqlon):
        """Find peak ground motion records associated with input event info.

        Ground motion records associated with input event are deleted from the
        database. Note that in the case of duplicate stations, the amps from
        only one will
        be used, any others will be deleted from the database.

        Args:
            eqtime (int): Unix timestamp of earthquake origin.
            eqlat (float): Latitude of earthquake origin.
            eqlon (float): Longitude of earthquake origin.
        Returns:
            list: A list of amps associated with the event. Each row in the
            list has the following columns:

                       - code: Station code
                       - channel: Channel (HHE, HHN, etc.)
                       - imt: Intensity measure type (pga, pgv, etc.)
                       - value: IMT value.
                       - lat: Station latitude.
                       - lon: Station longitude.
                       - netid: Station contributing network.
                       - name: String describing station name.
                       - distance: Distance (km) from station to origin.
                       - flag: Value will be 0.
                       - loccode: The location code of the instrument.
        """
        self._cursor.execute('BEGIN EXCLUSIVE')
        time_query = ('SELECT id, timestamp, lat, lon, code, network '
                      'FROM station WHERE timestamp > ? AND timestamp < ? ')
        self._cursor.execute(time_query, ((eqtime - TMIN), (eqtime + TMAX)))
        # numpy array of id, timestamp, lat, lon
        eqdata = np.array(self._cursor.fetchall())

        if not len(eqdata):
            self.commit()
            return []

        dist = geodetic_distance(eqlon, eqlat, eqdata[:, 3].astype(float),
                                 eqdata[:, 2].astype(float))
        inear = np.where(dist < DISTANCE)[0]
        eqdata = eqdata[inear]
        dist = dist[inear]
        stadict = {}
        junk_sids = []
        for idx, row in enumerate(eqdata):
            sid, timestamp, code, network = [row[x] for x in (0, 1, 4, 5)]
            timestamp = int(timestamp)
            if network not in stadict:
                stadict[network] = {code: {'sid': sid,
                                           'timestamp': timestamp,
                                           'distance': dist[idx]}}
                continue
            elif code not in stadict[network]:
                stadict[network][code] = {'sid': sid,
                                          'timestamp': timestamp,
                                          'distance': dist[idx]}
                continue
            traveltime = dist[idx] / 4.2
            new_dt = abs(abs(eqtime - timestamp) - traveltime)
            old_dt = abs(abs(eqtime - stadict[network][code]['timestamp'])
                         - traveltime)
            if old_dt < new_dt:
                junk_sids.append(sid)
                continue
            junk_sids.append(stadict[network][code]['sid'])
            stadict[network][code] = {'sid': sid,
                                      'timestamp': timestamp,
                                      'distance': dist[idx]}

        sta_sids = []
        for netd in stadict.values():
            for coded in netd.values():
                sta_sids.append(coded['sid'])

        if not len(sta_sids):
            self.commit()
            return []

        amp_query = ('SELECT s.network, s.name, s.code, s.lat, s.lon, '
                     'c.channel, c.loc, p.imt, p.value FROM station s, '
                     'channel c, pgm p WHERE s.id IN %s AND '
                     'c.station_id = s.id AND p.channel_id = c.id '
                     'ORDER BY s.network, s.code, c.channel, p.imt')
        delete_query = 'DELETE FROM station where id in %s'

        # data_list will hold the rows of the dataframe
        nstas = len(sta_sids)
        data_list = []
        start = 0
        while start < nstas:
            end = start + MAX_VARS
            if end > nstas:
                end = nstas
            varstr = '({0})'.format(
                ', '.join('?' for _ in sta_sids[start:end]))
            query = amp_query % varstr
            self._cursor.execute(query, sta_sids[start:end])
            amprows = self._cursor.fetchall()
            for row in amprows:
                # data_row = (code, channel_name, imt, value, lat, lon,
                #             network, name, distance, flag, loccode)
                data_row = (row[2], row[5], row[7], row[8], row[3], row[4],
                            row[0], row[1],
                            stadict[row[0]][row[2]]['distance'], 0, row[6])
                data_list.append(data_row)
            # Delete the stations now, since we have them queued up
            self._cursor.execute(delete_query % varstr, sta_sids[start:end])
            start = end

        # clean up rows that have been associated but didn't make the cut
        start = 0
        njunk = len(junk_sids)
        while start < njunk:
            end = start + MAX_VARS
            if end > njunk:
                end = njunk
            varstr = '({0})'.format(
                ', '.join('?' for _ in junk_sids[start:end]))
            self._cursor.execute(delete_query % varstr, junk_sids[start:end])
            start = end

        self.commit()
        return data_list

    def writeXML(self, data_list, eventid, pretty_print=False):
        """Write the list of tuples as an XML file in the event's
        current directory.

        Args:
            data_list (list): A list of tuples with the following
                              elements:

                              - station code
                              - channel
                              - imt
                              - imt value
                              - station latitude
                              - station longitude
                              - station's network id
                              - station's name string
                              - distance from station to origin
                              - imt flag
                              - channel's location code

            eventid (str): The event ID of the event associated with the data.
            pretty_print (bool): Whether or not to write the XML in a more
                                 human-readable form. If True, the file will
                                 be somewhat larger and writing will be
                                 somewhat slower.

        Returns:
            nothing: Nothing.
        """
        root = ET.Element('shakemap-data', code_version="4.0")
        create_time = int(time.time())
        stationlist = ET.SubElement(root, 'stationlist',
                                    created='%i' % create_time)
        oldnet = None
        oldcode = None
        oldchan = None
        oldloc = None
        for row in data_list:
            code, chan, imt, value, lat, lon, net, name, dist, flag, loc = row
            if net != oldnet or code != oldcode:
                if not code.startswith(net + '.'):
                    stacode = net + '.' + code
                else:
                    stacode = code
                station = ET.SubElement(stationlist, 'station', code=stacode,
                                        name=name, insttype="",
                                        lat="%.4f" % lat,
                                        lon="%.4f" % lon,
                                        dist="%.4f" % dist, netid=net,
                                        commtype="DIG", loc="")
                oldnet = net
                oldcode = code
                oldchan = None
                oldloc = None
            if chan != oldchan or loc != oldloc:
                if not chan.startswith(loc + '.'):
                    comp = loc + '.' + chan
                else:
                    comp = chan
                component = ET.SubElement(station, 'comp', name=comp)
                oldchan = chan
                oldloc = loc
            ET.SubElement(component, imt, value="%.6f" % value,
                          flag="%s" % str(flag))

        data_folder = os.path.join(self._data_path, eventid, 'current')
        if not os.path.isdir(data_folder):
            os.makedirs(data_folder)
        amptime = datetime.utcnow().strftime('%Y%m%d%H%M%S')
        xmlfile = os.path.join(data_folder, 'unassoc_%s_dat.xml' % amptime)

        if pretty_print:
            pstring = prettify(root)
            with open(xmlfile, 'w') as fd:
                fd.write(pstring)
        else:
            tree = ET.ElementTree(root)
            tree.write(xmlfile, encoding='utf-8', xml_declaration=True)

        return

    def __del__(self):
        """Destructor.

        """
        if hasattr(self, '_connection') and self._connection is not None:
            self._disconnect()

    def insertAmps(self, xmlfile):
        """Insert data from amps file into database.

        Args:
            xmlfile (str): XML file containing peak ground motion data.
        """
        _, fname = os.path.split(xmlfile)
        try:
            xmlstr = open(xmlfile, 'r').read()
            # sometimes these records have non-ascii bytes in them
            newxmlstr = re.sub(r'[^\x00-\x7F]+', ' ', xmlstr)
            # newxmlstr = _invalid_xml_remove(xmlstr)
            newxmlstr = newxmlstr.encode('utf-8', errors='xmlcharrefreplace')
            amps = dET.fromstring(newxmlstr)
        except Exception as e:
            raise Exception('Could not parse %s, due to error "%s"' %
                            (xmlfile, str(e)))

        if amps.tag != 'amplitudes':
            raise Exception('%s does not appear to be an amplitude XML '
                            'file.' % xmlfile)
        agency = amps.get('agency')
        record = amps.find('record')
        timing = record.find('timing')
        reference = timing.find('reference')
        has_pgm = False
        time_dict = {}
        for child in reference.iter():
            node_name = child.tag
            if node_name == 'PGMTime':
                has_pgm = True
            elif node_name == 'year':
                time_dict['year'] = int(child.get('value'))
            elif node_name == 'month':
                time_dict['month'] = int(child.get('value'))
            elif node_name == 'day':
                time_dict['day'] = int(child.get('value'))
            elif node_name == 'hour':
                time_dict['hour'] = int(child.get('value'))
            elif node_name == 'minute':
                time_dict['minute'] = int(child.get('value'))
            elif node_name == 'second':
                time_dict['second'] = int(child.get('value'))
            elif node_name == 'msec':
                time_dict['msec'] = int(child.get('value'))
        if has_pgm:
            pgmtime_str = reference.find('PGMTime').text
            try:
                tfmt = constants.TIMEFMT.replace('Z', '')
                pgmdate = datetime.strptime(
                    pgmtime_str[0:19], tfmt).replace(tzinfo=timezone.utc)
            except ValueError:
                tfmt = constants.ALT_TIMEFMT.replace('Z', '')
                pgmdate = datetime.strptime(
                    pgmtime_str[0:19], tfmt).replace(tzinfo=timezone.utc)
            pgmtime = int(dt_to_timestamp(pgmdate))
        else:
            if not len(time_dict):
                print('No time data for file %s' % fname)
                return
            pgmdate = datetime(time_dict['year'],
                               time_dict['month'],
                               time_dict['day'],
                               time_dict['hour'],
                               time_dict['minute'],
                               time_dict['second'])
            pgmtime = dt_to_timestamp(pgmdate)

        # there are often multiple stations per file, but they're
        # all duplicates of each other, so just grab the information
        # from the first one
        station = record.find('station')
        attrib = dict(station.items())
        lat = float(attrib['lat'])
        lon = float(attrib['lon'])
        code = attrib['code']
        name = attrib['name']
        if 'net' in attrib:
            network = attrib['net']
        elif 'netid' in attrib:
            network = attrib['netid']
        else:
            network = agency
        #
        # The station (at this pgmtime +/- 10 seconds) might already exist
        # in the DB; if it does, use it
        #
        self._cursor.execute('BEGIN EXCLUSIVE')
        query = ('SELECT id, timestamp FROM station where network = ? and '
                 'code = ? and timestamp > ? and timestamp < ?')
        self._cursor.execute(query,
                             (network, code, pgmtime - 10, pgmtime + 10))
        #
        # It's possible that the query returned more than one station; pick
        # the one closest to the new station's pgmtime
        #
        rows = self._cursor.fetchall()
        best_sid = None
        best_time = None
        for row in rows:
            dtime = abs(row[1] - pgmtime)
            if best_time is None or dtime < best_time:
                best_time = dtime
                best_sid = row[0]
        inserted_station = False
        if best_sid is None:
            fmt = ('INSERT INTO station '
                   '(timestamp, lat, lon, name, code, network) '
                   'VALUES (?, ?, ?, ?, ?, ?)')
            self._cursor.execute(fmt, (pgmtime, lat, lon, name, code, network))
            best_sid = self._cursor.lastrowid
            inserted_station = True

        #
        # If the station is already there, it has at least one channel, too
        #
        existing_channels = {}
        if inserted_station is False:
            chan_query = 'SELECT channel, id FROM channel where station_id = ?'
            self._cursor.execute(chan_query, (best_sid,))
            rows = self._cursor.fetchall()
            existing_channels = dict(rows)

        # might need these
        insert_channel = ('INSERT INTO channel '
                          '(station_id, channel, loc)'
                          'VALUES (?, ?, ?)')
        insert_pgm = ('INSERT INTO pgm '
                      '(channel_id, imt, value)'
                      'VALUES (?, ?, ?)')

        # loop over components
        channels_inserted = 0
        for channel in record.iter('component'):
            # We don't want channels with qual > 4 (assuming qual is Cosmos
            # table 6 value)
            qual = channel.get('qual')
            if qual:
                try:
                    iqual = int(qual)
                except ValueError:
                    # qual is something we don't understand
                    iqual = 0
            else:
                iqual = 0
            if iqual > 4:
                continue
            loc = channel.get('loc')
            if not loc:
                loc = '--'
            cname = channel.get('name')
            if cname in existing_channels:
                best_cid = existing_channels[cname]
                inserted_channel = False
            else:
                self._cursor.execute(insert_channel, (best_sid, cname, loc))
                best_cid = self._cursor.lastrowid
                inserted_channel = True
                channels_inserted += 1

            #
            # Similarly, if the channel is already there, we don't want to
            # insert repeated IMTs (and updating them doesn't make a lot of
            # sense)
            #
            existing_pgms = {}
            if inserted_channel is False:
                pgm_query = 'SELECT imt, id FROM pgm where channel_id = ?'
                self._cursor.execute(pgm_query, (best_cid,))
                rows = self._cursor.fetchall()
                existing_pgms = dict(rows)
            # loop over imts in channel
            pgm_list = []
            for pgm in list(channel):
                imt = pgm.tag
                if imt not in IMTS:
                    continue
                try:
                    value = float(pgm.get('value'))
                except ValueError:
                    #
                    # Couldn't interpret the value for some reason
                    #
                    continue
                if imt == 'sa':
                    imt = 'p'+imt+pgm.get('period').replace('.', '')
                    value = value / 9.81
                if imt in IMTDICT:
                    imt = IMTDICT[imt]
                if imt == 'pga':
                    value = value / 9.81
                if imt in existing_pgms:
                    continue
                pgm_list.append((best_cid, imt, value))
            if len(pgm_list) > 0:
                #
                # Insert the new amps
                #
                self._cursor.executemany(insert_pgm, pgm_list)
            elif inserted_channel:
                #
                # If we didn't insert any amps, but we inserted the channel,
                # delete the channel
                #
                channel_delete = 'DELETE FROM channel WHERE id = ?'
                self._cursor.execute(channel_delete, (best_cid,))
                channels_inserted -= 1
            # End of pgm loop
        # End of channel loop

        #
        # If we inserted the station but no channels, delete the station
        #
        if channels_inserted == 0 and inserted_station:
            station_delete = 'DELETE FROM station WHERE id = ?'
            self._cursor.execute(station_delete, (best_sid,))

        self.commit()
        return

    def cleanAmps(self, threshold=30):
        """Clean out amplitude data that is older than the threshold
        number of days.

        Args:
            threshold (int): Maximum age in days of amplitude data in
                             the database.
        Returns:
            int: Number of stations deleted.
        """
        thresh_date = dt_to_timestamp(datetime.utcnow() -
                                      timedelta(days=threshold))
        squery = 'DELETE FROM station WHERE timestamp < ?'
        self._cursor.execute(squery, [thresh_date])
        nrows = self._cursor.rowcount
        self.commit()
        return nrows

    def cleanEvents(self, threshold=365):
        """Clean out event data that is older than the threshold number
        of days.

        Args:
            threshold (int): Maximum age in days of events in the database.
        Returns:
            int: Number of events deleted.
        """
        thresh_date = dt_to_timestamp(datetime.utcnow() -
                                      timedelta(days=threshold))
        equery = 'DELETE FROM event WHERE time < %i' % thresh_date
        self._cursor.execute(equery)
        nevents = self._cursor.rowcount
        self.commit()
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
                  - station_min: Datetime of earliest amplitude data in
                                 database.
                  - station_max: Datetime of most recent amplitude data
                                 in database.
        """
        results = {}

        # event stuff
        equery = 'SELECT count(*), min(time), max(time) FROM event'
        self._cursor.execute(equery)
        row = self._cursor.fetchone()
        results['events'] = row[0]
        if row[0] == 0:
            results['event_min'] = None
            results['event_max'] = None
        else:
            results['event_min'] = datetime.fromtimestamp(row[1],
                                                          timezone.utc)
            results['event_max'] = datetime.fromtimestamp(row[2],
                                                          timezone.utc)

        # station stuff
        squery = 'SELECT count(*), min(timestamp), max(timestamp) FROM station'
        self._cursor.execute(squery)
        row = self._cursor.fetchone()
        results['stations'] = row[0]
        if row[0] == 0:
            results['station_min'] = None
            results['station_max'] = None
        else:
            results['station_min'] = datetime.fromtimestamp(row[1],
                                                            timezone.utc)
            results['station_max'] = datetime.fromtimestamp(row[2],
                                                            timezone.utc)

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

        return results


def dt_to_timestamp(dt):
    timestamp = int(dt.replace(tzinfo=timezone.utc).timestamp())
    return timestamp


def timestr_to_timestamp(timestr):
    try:
        timestamp = int(datetime.strptime(timestr, constants.TIMEFMT).
                        replace(tzinfo=timezone.utc).timestamp())
    except ValueError:
        timestamp = int(datetime.strptime(timestr, constants.ALT_TIMEFMT).
                        replace(tzinfo=timezone.utc).timestamp())
    return timestamp


# def _invalid_xml_remove(c):
    # http://stackoverflow.com/questions/1707890/fast-way-to-filter-illegal-xml-unicode-chars-in-python
    # noqa
#    illegal_unichrs = [(0x00, 0x08), (0x0B, 0x1F), (0x7F, 0x84), (0x86, 0x9F),
#                       (0xD800, 0xDFFF), (0xFDD0, 0xFDDF), (0xFFFE, 0xFFFF),
#                       (0x1FFFE, 0x1FFFF), (0x2FFFE, 0x2FFFF),
#                       (0x3FFFE, 0x3FFFF), (0x4FFFE, 0x4FFFF),
#                       (0x5FFFE, 0x5FFFF), (0x6FFFE, 0x6FFFF),
#                       (0x7FFFE, 0x7FFFF), (0x8FFFE, 0x8FFFF),
#                       (0x9FFFE, 0x9FFFF), (0xAFFFE, 0xAFFFF),
#                       (0xBFFFE, 0xBFFFF), (0xCFFFE, 0xCFFFF),
#                       (0xDFFFE, 0xDFFFF), (0xEFFFE, 0xEFFFF),
#                       (0xFFFFE, 0xFFFFF),
#                       (0x10FFFE, 0x10FFFF)]
#
#    illegal_ranges = ["%s-%s" % (chr(low), chr(high))
#                      for (low, high) in illegal_unichrs
#                      if low < sys.maxunicode]
#
#    illegal_xml_re = re.compile(u'[%s]' % u''.join(illegal_ranges))
#    if illegal_xml_re.search(c) is not None:
#        # Replace with space
#        return ' '
#    else:
#        return c


def prettify(elem):
    """Return a pretty-printed XML string.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
