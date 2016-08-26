#!/usr/bin/env python

# stdlib imports
import sqlite3
from xml.dom import minidom
import os.path
import sys
import copy
import time
from collections import OrderedDict

# third party imports
import pandas as pd
import numpy as np
from openquake.hazardlib import imt as GEM_IMT
from openquake.hazardlib.gsim import base
from openquake.hazardlib import const

# local imports
from shakemap.grind.gmice.wgrw12 import WGRW12
from .distance import get_distance, get_distance_measures

TABLES = {'station':
          {'id': 'integer primary key',
           'network': 'str',
           'code': 'str',
           'name': 'str',
           'lat': 'float',
           'lon': 'float',
           'elev': 'float',
           'vs30': 'float',
           'instrumented': 'int'},  # distance colums will be added when the table is created
          'imt':
          {'id': 'integer primary key',
           'imt_type': 'str'},
          'amp':
          {'id': 'integer primary key',
           'station_id': 'int',
           'imt_id': 'int',
           'original_channel': 'str',
           'orientation': 'str',
           'amp': 'float',
           'uncertainty': 'float',
           'flag': 'str'},
          'soiltype':
          {'id': 'integer primary key',
           'soil_type': 'str'},
          'predicted':
          {'id': 'integer primary key',
           'station_id': 'int',
           'imt_id': 'int',
           'soiltype_id': 'int',
           'amp': 'float',
           'uncertainty_total': 'float',
           'uncertainty_inter': 'float',
           'uncertainty_intra': 'float'},
          'siteamp':
          {'id': 'integer primary key',
           'station_id': 'int',
           'imt_id': 'int',
           'amp_factor': 'float'}
          }

BASE_IMTS = ['mmi',
             'pga',
             'pgv',
             'psa03',
             'psa10',
             'psa30']

IMT_TYPES = {'mmi': 0, 
             'pga': 1, 
             'pgv': 2, 
             'psa03': 3, 
             'psa10': 4, 
             'psa30': 5,
             'mmi_from_pga': 6, 
             'mmi_from_pgv': 7, 
             'mmi_from_psa03': 8, 
             'mmi_from_psa10': 9, 
             'mmi_from_psa30': 10,
             'pga_from_mmi': 11, 
             'pgv_from_mmi': 12, 
             'psa03_from_mmi': 13, 
             'psa10_from_mmi': 14, 
             'psa30_from_mmi': 15 }

IMT_TYPES_ORDERED = OrderedDict(sorted(IMT_TYPES.items(), key=lambda t: t[1]))

IMT_LOOKUP = { value: key for key, value in IMT_TYPES.items() }

# dictionary of our imt type strings to the kind that GEM needs to create
# IMT objects.
GEM_IMT_MAP = {'mmi': 'MMI',
               'pga': 'PGA',
               'pgv': 'PGV',
               'psa03': 'SA(0.3)',
               'psa10': 'SA(1.0)',
               'psa30': 'SA(3.0)'}

SOIL_TYPES = {'rock': 0, 'soil': 1}

DISTANCES = get_distance_measures()

def _getOrientation(orig_channel):
    if orig_channel[-1] in ('N', 'E', 'Z'):
        orientation = orig_channel[-1]
    else:
        orientation = 'U'  # this is unknown
    return orientation


def _getGroundMotions(comp):
    """
    Get a dictionary of peak ground motions (values and flags).  Output keys are one of: [pga,pgv,psa03,psa10,psa30]
    Even if flags are not specified in the input, they will be guaranteed to at least have a flag of '0'.
    """
    pgmdict = {}
    for pgm in comp.childNodes:
        if pgm.nodeName == '#text':
            continue
        key = pgm.nodeName
        if key == 'acc':
            key = 'pga'
        if key == 'vel':
            key = 'pgv'

        value = float(pgm.getAttribute('value'))
        if np.isnan(value):
            x = 1
        if pgm.hasAttribute('flag'):
            flag = pgm.getAttribute('flag')
        else:
            flag = '0'
        pgmdict[key] = {'value': value, 'flag': flag}
    return pgmdict


def _getStationAttributes(station):
    """
    Get a dictionary of the station attributes
    """
    attrdict = {}
    for attr in list(station.attributes.items()):
        key = attr[0]
        value = attr[1]
        # is this value a float or str?
        try:
            value = float(value)
        except:
            pass
        attrdict[key] = value
    return attrdict


def _filter_station(xmlfile):
    """
    Filter individual xmlfile into a stationdict data structure.
    Inputs:

     * xmlfile xml file (or file-like object) containing station data

    Outputs:

     * stationdict Data structure as returned by filter_stations()
    """
    stationdict = {}
    dom = minidom.parse(xmlfile)
    for root in dom.childNodes:
        if not isinstance(root, minidom.DocumentType):
            break
    stations = root.getElementsByTagName('station')
    for station in stations:
        code = station.getAttribute('code')
        attributes = _getStationAttributes(station)
        comps = station.getElementsByTagName('comp')
        if code in list(stationdict.keys()):
            compdict = stationdict[code]
        else:
            compdict = {}
        for comp in comps:
            compname = comp.getAttribute('name')
            tpgmdict = _getGroundMotions(comp)
            if compname in list(compdict.keys()):
                pgmdict = compdict[compname]
            else:
                pgmdict = {}
            pgmdict.update(tpgmdict)
            # copy the VALUES, not REFERENCES, of the component list into our
            # growing dictionary
            compdict[compname] = copy.deepcopy(pgmdict)
        if 'intensity' in list(attributes.keys()):
            compdict[compname]['mmi'] = {
                'value': attributes['intensity'], 'flag': '0'}
        stationdict[code] = (attributes, copy.deepcopy(compdict))
    dom.unlink()
    return stationdict


def _createTables(db, cursor):
    for table in TABLES.keys():
        sql = 'CREATE TABLE %s (' % table
        nuggets = []
        for column, ctype in TABLES[table].items():
            nuggets.append('%s %s' % (column, ctype))
        sql += ','.join(nuggets) + ')'
        cursor.execute(sql)
        db.commit()

    # Add distance measures to the station table
    for col in DISTANCES:
        cursor.execute('ALTER TABLE station ADD COLUMN %s float' % (col))
    db.commit()

    # IMT types are either observed (first row here)
    # derived MMI (second row)
    # or derived PGM (third row)
    rows = []
    for imt_type, imt_id in IMT_TYPES.items():
        rows.append((imt_id, imt_type))
    cursor.executemany('INSERT INTO imt (id, imt_type) VALUES (?, ?)', rows)
    db.commit()
    for soiltype, sid in SOIL_TYPES.items():
        soilquery = 'INSERT INTO soiltype (soil_type, id) VALUES ("%s", %d)' % (soiltype, sid)
        cursor.execute(soilquery)
        db.commit()


class StationList(object):

    CIIM_TUPLE = ('dyfi', 'mmi', 'intensity', 'ciim')

    def __init__(self, dbfile):
        self.db = sqlite3.connect(dbfile)
        self.cursor = self.db.cursor()


    def __len__(self):
        squery = 'SELECT count(*) FROM station'
        self.cursor.execute(squery)
        return self.cursor.fetchone()[0]


    def __del__(self):
        self.cursor.close()
        self.db.close()


    @classmethod
    def loadFromXML(cls, xmlfiles, dbfile):
        stationdictlist = []
        for xmlfile in xmlfiles:
            stationdict = _filter_station(xmlfile)
            stationdictlist.append(stationdict)
        return cls.loadFromDict(stationdictlist, dbfile)


    @classmethod
    def loadFromDict(cls, stationdictlist, dbfile):
        do_create = False
        if not os.path.isfile(dbfile):
            do_create = True
        dbfile = dbfile
        db = sqlite3.connect(dbfile)
        cursor = db.cursor()
        if do_create:
            # create the tables we want
            _createTables(db, cursor)

        sid = 0
        amp_rows = []
        station_rows = []
        for stationdict in stationdictlist:
            for key, station_tpl in stationdict.items():
                station_attributes, comp_dict = station_tpl
                lat = station_attributes['lat']
                lon = station_attributes['lon']
                network = station_attributes['netid']
                code = key
                if key.startswith(network):
                    code = key.replace(network + '.', '')
                name = station_attributes['name']
                # elevation?

                #
                # How to determine a station is instrumented?
                #
                instrumented = int(station_attributes['netid'].lower() 
                                   not in cls.CIIM_TUPLE)
                station_rows.append((sid, network, code, name, lat, lon, 
                                     instrumented))

                for original_channel, pgm_dict in comp_dict.items():
                    orientation = _getOrientation(original_channel)
                    for imt_type, imt_dict in pgm_dict.items():
                        if imt_type not in IMT_TYPES:
                            continue
                        if (instrumented == 0) and (imt_type != 'mmi'):
                            continue
                        imtid = IMT_TYPES[imt_type]
                        amp = imt_dict['value']
                        flag = imt_dict['flag']

                        if np.isnan(amp):
                            amp = 'NULL'
                        amp_rows.append((sid, imtid, original_channel, 
                                         orientation, amp, flag))
                sid += 1

        fmt = 'INSERT INTO station (id, network, code, name, lat, lon, '\
              'instrumented) VALUES (?, ?, ?, ?, ?, ?, ?)'
        cursor.executemany(fmt, station_rows)
        db.commit()

        fmt = 'INSERT INTO amp (station_id, imt_id, original_channel, '\
              'orientation, amp, flag) VALUES (?, ?, ?, ?, ?, ?)'
        cursor.executemany(fmt, amp_rows)
        db.commit()
        cursor.close()
        db.close()
        return cls(dbfile)


    def fillTables(self, source, vs30, gmpe, gmice=None, ipe=None):
        """
        Populate tables with derived MMI/PGM values and distances.

        :param source:
          ShakeMap Source object.
        """
        gmice = WGRW12()
        emag = source.getEventDict()['mag']
        #
        # Do the distances for all of the stations
        #
        query = 'SELECT id, lat, lon, code, network FROM station'
        self.cursor.execute(query)
        station_rows = self.cursor.fetchall()

        nrows = len(station_rows)
        lats = np.empty((nrows))
        lons = np.empty((nrows))
        depths = np.zeros((nrows))
        for irow in range(nrows):
            lats[irow] = station_rows[irow][1]
            lons[irow] = station_rows[irow][2]
        ddict = get_distance(DISTANCES, lats, lons, depths, source)
        dist_rows = []
        for irow in range(nrows):
            dist_rows.append(tuple(ddict[dt][irow] for dt in DISTANCES) + (station_rows[irow][0],))

        query = 'UPDATE station set '
        for dt in DISTANCES:
            query += dt + '=?'
            if dt == DISTANCES[-1]:
                query += ' '
            else:
                query += ', '
        query += 'WHERE id=?'
        self.cursor.executemany(query, dist_rows)
        self.db.commit()

        # 
        # Use the GMICE to get the estimated MMIs from the various 
        # ground motion parameters from the instrumented stations.
        # These values will fill the rows like 'mmi_from_pgm'
        #
        amp_rows = []
        for imt, imtid in IMT_TYPES_ORDERED.items():
            if imt.endswith('_mmi') or imt.startswith('mmi'):
                continue
            query = 'SELECT a.original_channel, a.orientation, a.uncertainty, '\
                    'a.amp, s.repi, s.id FROM amp a, station s '\
                    'WHERE a.imt_id=%i AND a.station_id=s.id AND '\
                    's.instrumented=1 AND a.flag=0' % (imtid)
            self.cursor.execute(query)
            rows = self.cursor.fetchall()

            nrows = len(rows)
            amps = np.empty((nrows))
            dists = np.empty((nrows))
            emags = np.zeros((nrows)) + emag
            for irow in range(nrows):
                amps[irow] = rows[irow][3]
                dists[irow] = rows[irow][4]

            gemimt = GEM_IMT.from_string(GEM_IMT_MAP[imt])
            dmmi = gmice.getMIfromGM(amps, gemimt, dists=dists, mag=emags)

            derived_imtid = IMT_TYPES["mmi_from_" + imt]

            for irow in range(nrows):
                amp_rows.append((rows[irow][5], derived_imtid, rows[irow][0],
                                 rows[irow][1], dmmi[irow], 0.0))

        self.cursor.executemany(
            'INSERT INTO amp (station_id, imt_id, original_channel, '\
            'orientation, amp, uncertainty, flag) VALUES '\
            '(?, ?, ?, ?, ?, ?, "0")', amp_rows)

        # 
        # Use the GMICE to get the estimated ground motions from the 
        # mmi from the non-instrumented stations.
        # These will appeara in columns named like 'pgm_from_mmi'
        #
        imtid = IMT_TYPES['mmi']
        query = 'SELECT a.uncertainty, a.amp, s.repi, s.id FROM amp a, '\
                'station s WHERE a.imt_id=%i AND a.station_id=s.id AND '\
                's.instrumented=0 AND a.flag=0' % (imtid)
        self.cursor.execute(query)
        rows = self.cursor.fetchall()

        nrows = len(rows)
        mmi = np.empty((nrows))
        dists = np.empty((nrows))
        emags = np.zeros((nrows)) + emag
        for irow in range(nrows):
            mmi[irow] = rows[irow][1]
            dists[irow] = rows[irow][2]

        amp_rows = []
        for imt, imtid in IMT_TYPES_ORDERED.items():
            if imt.endswith('_mmi') or imt.startswith('mmi'):
                continue

            gemimt = GEM_IMT.from_string(GEM_IMT_MAP[imt])
            dmmi = gmice.getGMfromMI(mmi, gemimt, dists=dists, mag=emags)

            derived_imtid = IMT_TYPES[imt + "_from_mmi"]

            for irow in range(nrows):
                amp_rows.append((rows[irow][3], derived_imtid, dmmi[irow], 
                                 0.0))

        self.cursor.executemany(
            'INSERT INTO amp (station_id, imt_id, amp, uncertainty, flag) '\
            'VALUES (?, ?, ?, ?, "0")', amp_rows)
        self.db.commit()

        #
        # Use the gmpe to predict ground motions at each station
        #
        rx = source.getRuptureContext([gmpe])
        dx = base.DistancesContext()
        for method in DISTANCES:
            (dx.__dict__)[method] = ddict[method]
        sx = base.SitesContext()
        sx.lats = lats
        sx.lons = lons
        vs30_soil = vs30.getValue(lats, lons, default=760.0)
        vs30_rock = np.zeros_like(vs30_soil) + 760.0
#        sx.z1pt0 = sites._calculate_z1pt0(sx.vs30)
#        sx.z2pt5 = sites._calculate_z2pt5(sx.z1pt0)
        stddev_types = [const.StdDev.TOTAL, const.StdDev.INTER_EVENT, 
                        const.StdDev.INTRA_EVENT]
        pred_rows = []
        siteamp_rows = []
        vs30_rows = []
        nrows = len(station_rows)
        for imt in BASE_IMTS:
            gemimt = GEM_IMT.from_string(GEM_IMT_MAP[imt])
            sx.vs30 = vs30_soil
            pred_soil, pred_stdev = gmpe.get_mean_and_stddevs(sx, rx, dx, imt, stddev_types)
            sx.vs30 = vs30_rock
            pred_rock, junk = gmpe.get_mean_and_stddevs(sx, rx, dx, imt, stddev_types)
            amp_facts = pred_soil - pred_rock

            for irow in range(nrows):
                pred_rows.append(
                    (station_rows[irow][0], IMT_TYPES[imt], SOIL_TYPES['rock'], 
                     pred_rock[irow], pred_stdev[const.StdDev.TOTAL][irow],
                     pred_stdev[const.StdDev.INTER_EVENT][irow],
                     pred_stdev[const.StdDev.INTRA_EVENT][irow]))
                siteamp_rows.append(
                    (station_rows[irow][0], IMT_TYPES[imt], amp_facts[irow]))

        for irow in range(nrows):
            vs30_rows.append((vs30_soil[irow], station_rows[irow][0]))
            
        self.cursor.executemany(
            'INSERT INTO predicted (station_id, imt_id, soiltype_id, amp, '\
            'uncertainty_total, uncertainty_inter, uncertainty_intra) '\
            'VALUES (?, ?, ?, ?, ?, ?, ?)', pred_rows)
        self.cursor.executemany(
            'INSERT INTO siteamp (station_id, imt_id, amp_factor) '\
            'VALUES (?, ?, ?)', siteamp_rows)
        self.cursor.executemany(
            'UPDATE station SET vs30=? WHERE id=?', vs30_rows)

        self.db.commit()


    def getInstrumentedStations(self):
        dstr = ''
        columns = ['id', 'lat', 'lon', 'code', 'network']
        for mm in DISTANCES:
            dstr += ", %s" % (mm)
            columns.append(mm)

        stationquery = 'SELECT id, lat, lon, code, network%s FROM station '\
                       'where instrumented = 1' % (dstr)
        
        self.cursor.execute(stationquery)
        rows = self.cursor.fetchall()
        nrows = len(rows)
        try:
            df = pd.DataFrame(rows, columns=columns)
        except e:
            print('Exception in creating data frame')
            raise e

        ncols = 0
        new_cols = []
        for imt in IMT_TYPES_ORDERED.keys():
            new_cols.append(imt)
            new_cols.append(imt + '_unc')
            ncols += 2
        app = np.empty((np.shape(rows)[0], ncols))
        app[:] = np.nan

        col_dict = dict(zip(new_cols, range(ncols)))
        id_dict = dict(zip(df['id'], range(nrows)))

        #
        # Get all of the unflagged instrumented amps with the proper
        # orientation
        #
        self.cursor.execute(
            'SELECT a.amp, a.uncertainty, a.imt_id, a.station_id FROM '\
            'amp a, station s WHERE a.flag = "0" AND s.id = a.station_id '\
            'AND s.instrumented = 1 AND a.orientation NOT IN ("Z", "U") '\
            'AND a.amp IS NOT NULL'
            )
        amp_rows = self.cursor.fetchall()

        #
        # Go through all the instrumented amps and put them into the data 
        # frame
        #
        for this_row in amp_rows:
            imt_id = this_row[2]
            if imt_id not in IMT_LOOKUP:
                continue
            #
            # Set the cell to the peak amp
            #
            imt = IMT_LOOKUP[imt_id]
            rowidx = id_dict[this_row[3]]
            cval = app[rowidx, col_dict[imt]]
            amp = this_row[0]
            if np.isnan(cval) or (cval < amp):
                app[rowidx, col_dict[imt]] = amp
                app[rowidx, col_dict[imt + '_unc']] = this_row[1]

        df2 = pd.DataFrame(app, columns=new_cols)
        df = pd.concat([df, df2], axis=1)

        df['name'] = df.network.map(str) + '.' + df.code.map(str)
        del df['network']
        del df['code']
        newcols = ['name', 'lat', 'lon', ] + DISTANCES
        for imt in IMT_TYPES_ORDERED.keys():
            if not imt.startswith('mmi_from_'):
                continue
            newcols.append(imt)
            newcols.append(imt + '_unc')
        if pd.__version__ >= '0.17.0':
            df = df[newcols].sort_values('name')
        else:
            df = df[newcols].sort('name')
        return df


    def getMMIStations(self):
        dstr = ''
        columns = ['id', 'lat', 'lon', 'code', 'network']
        for mm in DISTANCES:
            dstr += ", %s" % (mm)
            columns.append(mm)

        stationquery = 'SELECT id, lat, lon, code, network%s FROM station '\
                       'where instrumented = 0' % (dstr)
        
        self.cursor.execute(stationquery)
        rows = self.cursor.fetchall()
        nrows = len(rows)
        try:
            df = pd.DataFrame(rows, columns=columns)
        except e:
            print('Exception in creating data frame')
            raise e

        ncols = 0
        new_cols = []
        for imt in IMT_TYPES_ORDERED.keys():
            new_cols.append(imt)
            new_cols.append(imt + '_unc')
            ncols += 2
        app = np.empty((np.shape(rows)[0], ncols))
        app[:] = np.nan

        col_dict = dict(zip(new_cols, range(ncols)))
        id_dict = dict(zip(df['id'], range(nrows)))

        #
        # Get all of the unflagged uninstrumented amps 
        #
        self.cursor.execute(
            'SELECT a.amp, a.uncertainty, a.imt_id, a.station_id FROM '\
            'amp a, station s WHERE a.flag = "0" AND s.id = a.station_id '\
            'AND s.instrumented = 0 '\
            'AND a.amp IS NOT NULL'
            )
        amp_rows = self.cursor.fetchall()

        #
        # Go through all the instrumented amps and put them into the data 
        # frame
        #
        for this_row in amp_rows:
            imt_id = this_row[2]
            if imt_id not in IMT_LOOKUP:
                continue
            #
            # Set the cell to the peak amp
            #
            imt = IMT_LOOKUP[imt_id]
            rowidx = id_dict[this_row[3]]
            cval = app[rowidx, col_dict[imt]]
            amp = this_row[0]
            if np.isnan(cval) or (cval < amp):
                app[rowidx, col_dict[imt]] = amp
                app[rowidx, col_dict[imt + '_unc']] = this_row[1]

        df2 = pd.DataFrame(app, columns=new_cols)
        df = pd.concat([df, df2], axis=1)

        df['name'] = df.network.map(str) + '.' + df.code.map(str)
        del df['network']
        del df['code']
        newcols = ['name', 'lat', 'lon', ] + DISTANCES
        for imt in IMT_TYPES_ORDERED.keys():
            if not imt.endswith('_from_mmi'):
                continue
            newcols.append(imt)
            newcols.append(imt + '_unc')
        if pd.__version__ >= '0.17.0':
            df = df[newcols].sort_values('name')
        else:
            df = df[newcols].sort('name')
        return df


if __name__ == '__main__':
    xmlfiles = sys.argv[1:]
    dbfile = 'stations.db'
    if os.path.isfile(dbfile):
        os.remove(dbfile)
    t1 = time.time()
    stations = StationList.loadFromXML(dbfile, xmlfiles)
    t2 = time.time()
    print('%i stations loaded in %.2f seconds' % (len(stations), t2 - t1))
    mmidf = stations.getMMIStations()
    t3 = time.time()
    print('%i intensity observations retrieved in %.2f seconds' %
          (len(mmidf), t3 - t2))
    imtdf = stations.getInstrumentedStations()
    t4 = time.time()
    print('%i instrumental measurements retrieved in %.2f seconds' %
          (len(imtdf), t4 - t3))
