# stdlib imports
import sqlite3
from xml.dom import minidom
import os.path
import copy
from collections import OrderedDict
import re

# third party imports
import pandas as pd
import numpy as np

# local imports
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
           'instrumented': 'int'},
          'imt':
          {'id': 'integer primary key',
           'type': 'str'},
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

SOIL_TYPES = {'rock': 0, 'soil': 1}

#
# These are the netid's that indicate MMI data
#
CIIM_TUPLE = ('dyfi', 'mmi', 'intensity', 'ciim')

#
# This is the full list of distance measure that we can calculate in the
# distance module. We compute all of them for all of the stations just
# because we can.
#
DISTANCES = get_distance_measures()


class StationList(object):
    """
    A class to facilitate reading ShakeMap formatted XML fies of peak
    amplitudes and MMI, and
    produce tables of station data. Seismic stations are considered to
    be 'instrumented'; MMI data is not instrumented and is indicated
    in the ShakeMap XML with a ``netid`` attribute of "DYFI," "MMI,"
    "INTENSITY," or "CIIM."

    .. note::
      Typically the user will call the class method :meth:`fromXML`
      to create a :class:`StationList` object the first time
      a set of station files are processed. (Or, as an alternative,
      the user can call :meth:`loadFromXML` and :meth:`fillTables`
      sequentially.)
      This will create a database at the location specified by the
      ``dbfile`` parameter to :meth:`fromXML`. Subsequent programs
      can use the default constructor to simply load ``dbfile``.

    """

    def __init__(self, dbfile):
        """
        The default constructor reads a pre-built SQLite database of
        station data.

        Args:
            dbfile (str):
                A SQLite database file containing pre-processed
                station data.

        Returns:
            A :class:`StationList` object.

        """
        self.IMT_TYPES = {}

        self.db = sqlite3.connect(dbfile)
        self.cursor = self.db.cursor()

        #
        # Fill in the IMT type lists/dicts
        #
        query = 'SELECT id, type FROM imt'
        self.cursor.execute(query)
        imt_rows = self.cursor.fetchall()

        for row in imt_rows:
            self.IMT_TYPES[row[1]] = row[0]

        #
        # This flips the key/value pairs in the IMT_TYPES dictionary
        #
        self.IMT_LOOKUP = {value: key for key, value in self.IMT_TYPES.items()}

        #
        # This is an ordered version of IMT_TYPES just to help make the tables
        # a little more human readable, and to keep them consistent for testing.
        #
        self.IMT_TYPES_ORDERED = OrderedDict(sorted(self.IMT_TYPES.items(), key=imt_sort))


    def __len__(self):
        """
        Returns the number of stations in the database.
        """
        squery = 'SELECT count(*) FROM station'
        self.cursor.execute(squery)
        return self.cursor.fetchone()[0]


    def __del__(self):
        """
        Closes out the database when the object is destroyed.
        """
        self.db.commit()
        self.cursor.close()
        self.db.close()


    @classmethod
    def fromXML(cls, xmlfiles, dbfile, origin, sites, rupture):
        """
        Create a StationList object by reading one or more ShakeMap XML
        input files and populate database tables with derived MMI/PGM
        values and distances. This is a convenience method that calls
        :meth:`loadFromXML` and :meth:`fillTables`.

        Args:
            xmlfiles (sequence of strings):
                Sequence of ShakeMap XML input files to read.
            dbfile (string):
                Path to a file into which to write the SQLite database. If
                the file exists, it will first be deleted.
            origin:
                ShakeMap Origin object containing information about the
                origin and source of the earthquake.
            sites:
                ShakeMap Sites object containing grid of Vs30 values for the
                region in question.
            rupture:
                ShakeMap Rupture object.

        Returns:
            :class:`StationList` object
        """

        #
        # We're creating a new database from the XML inputs, so
        # delete the existing dbfile if it exsits (yes, I can imagine
        # use cases where you might want to add to an existing
        # database, but we don't support that because our indexing
        # scheme doesn't currently allow it.)
        #
        try:
            os.remove(dbfile)
        except OSError:
            pass

        this = cls.loadFromXML(xmlfiles, dbfile)
        this.fillTables(origin, sites, rupture)
        return this


    @classmethod
    def loadFromXML(cls, xmlfiles, dbfile):
        """
        Create a StationList object by reading one or more ShakeMap XML input
        files.

        Args:
            xmlfiles (sequence of str):
                Sequence of ShakeMap XML input files to read.
            dbfile (str):
                Path to a file into which to write the SQLite database.

        Returns:
            :class:`StationList` object

        """

        stationdictlist = []
        imtset = set()
        for xmlfile in xmlfiles:
            stationdict, ims = cls._filter_station(xmlfile)
            stationdictlist.append(stationdict)
            imtset |= ims
        this = cls._loadFromDict(stationdictlist, dbfile, imtset)
        this.imtset = imtset
        return this


    @classmethod
    def _loadFromDict(cls, stationdictlist, dbfile, imtset):
        """
        Internal method to turn the station dictionary created from the
        ShakeMap XML input files into a SQLite database.

        Args:
            stationdictlist (list of stationdicts):
                A list of station dictionaries returned by _filter_station().
            dbfile (string):
                The path to which the SQLite database will be written.

        Returns:
            :class:`StationList` object
        """
        db = sqlite3.connect(dbfile)
        cursor = db.cursor()
        # create the tables we want
        sorted_imtset = sorted(list(imtset), key=imt_sort)
        cls._createTables(db, cursor, sorted_imtset)

        query = 'SELECT id, type FROM imt'
        cursor.execute(query)
        imt_rows = cursor.fetchall()

        imt_types = {}
        for row in imt_rows:
            imt_types[row[1]] = row[0]

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
                                   not in CIIM_TUPLE)

                station_rows.append((sid, network, code, name, lat, lon,
                                     instrumented))

                for original_channel, pgm_dict in comp_dict.items():
                    orientation = cls._getOrientation(original_channel)
                    for imt_type, imt_dict in pgm_dict.items():
                        if imt_type not in imt_types:
                            continue
                        if (instrumented == 0) and (imt_type != 'MMI'):
                            continue
                        imtid = imt_types[imt_type]
                        amp = imt_dict['value']
                        flag = imt_dict['flag']
                        if np.isnan(amp) or (amp <= 0):
                            amp = 'NULL'
                            flag = 'G'
                        elif imt_type == 'MMI':
                            pass
                        elif imt_type == 'PGV':
                            amp = np.log(amp)
                        else:
                            amp = np.log(amp / 100.0)

                        amp_rows.append((sid, imtid, original_channel,
                                         orientation, amp, flag))
                sid += 1

        cursor.executemany(
                'INSERT INTO station (id, network, code, name, lat, lon, '
                'instrumented) VALUES (?, ?, ?, ?, ?, ?, ?)', station_rows
            )

        cursor.executemany(
                'INSERT INTO amp (station_id, imt_id, original_channel, '
                'orientation, amp, flag) VALUES (?, ?, ?, ?, ?, ?)', amp_rows
            )
        db.commit()
        cursor.close()
        db.close()
        return cls(dbfile)


    def fillTables(self, origin, sites, rupture):
        """
        Populate database tables with derived MMI/PGM values and distances.
        This method should be called after :meth:`loadFromXML`.

        Args:
            origin:
                ShakeMap Origin object containing information about the
                origin and source of the earthquake.
            sites:
                ShakeMap Sites object containing grid of Vs30 values for the
                region in question.
            rupture:
                ShakeMap Rupture object.

        Returns:
            nothing
        """
        emag = origin.mag
        #
        # Get a list of stations
        #
        query = 'SELECT id, lat, lon, code, network FROM station'
        self.cursor.execute(query)
        station_rows = self.cursor.fetchall()

        #
        # Do the distances for all of the stations
        #
        nrows = len(station_rows)
        lats = np.empty((nrows))
        lons = np.empty((nrows))
        depths = np.zeros((nrows))
        for irow, row in enumerate(station_rows):
            lats[irow] = row[1]
            lons[irow] = row[2]
        ddict = get_distance(DISTANCES, lats, lons, depths, rupture)
        dist_rows = []
        for irow, row in enumerate(station_rows):
            dist_rows.append(
                    tuple(ddict[dt][irow] for dt in DISTANCES) +
                    (row[0],)
                )

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
        # store the Vs30 for each station
        #
        lldict = {'lats': lats, 'lons': lons}
        sx_soil = sites.getSitesContext(lldict)

        vs30_rows = []
        for irow, row in enumerate(station_rows):
            vs30_rows.append((sx_soil.vs30[irow], row[0]))

        self.cursor.executemany(
                'UPDATE station SET vs30=? WHERE id=?', vs30_rows)

        self.db.commit()


    def getStationDataframe(self, instrumented, sort=False):
        """
        Return a Pandas dataframe of the instrumented or non-instrumented
        stations.

        For the standard set of ShakeMap IMTs (mmi, pga, pgv, psa03, psa10,
        psa30), the columns in the dataframe would be:

        'id', 'lat', 'lon', 'code', 'network', 'vs30', 'repi', 'rhypo', 'rjb',
        'rrup', 'rx', 'ry', 'ry0', 'U', 'T', 'SA(0.3)', 'SA(3.0)', 'PGA', 'PGV',
        'SA(1.0)', 'name'

        For the non-instrumented dataframe, the columns would be:

        'id', 'lat', 'lon', 'code', 'network', 'vs30', 'repi', 'rhypo', 'rjb',
        'rrup', 'rx', 'ry', 'ry0', 'U', 'T', 'MMI', 'name'

        The **name** column is **network** and **code** concatenated with a
        period (".") between them.
        All ground motion
        units are natural log units. Distances are in km.


        Args:
            instrumented (integer):
                Set to 1 (one) if the dataframe is to contain the instrumented
                stations, or to 0 (zero) if the dataframe is to contain the
                non-instrumented (MMI) stations.
            sort (bool):
                If True, the dataframe will be sorted by the **name** column.
                The default if False (unsorted).

        Returns:
            A Pandas dataframe.
        """

        dstr = ''
        columns = []
        basic_columns = ['id', 'lat', 'lon', 'code', 'network', 'vs30']
        for mm in basic_columns + DISTANCES:
            if mm == basic_columns[0]:
                dstr = '%s' % (mm)
            else:
                dstr += ", %s" % (mm)
            columns.append(mm)

        self.cursor.execute(
                'SELECT %s FROM station where instrumented = %d' %
                (dstr, instrumented)
            )

        station_rows = self.cursor.fetchall()
        nstation_rows = len(station_rows)
        try:
            df = pd.DataFrame(station_rows, columns=columns)
        except:
            print('Exception in creating data frame')
            raise

        ncols_amp = 0
        new_cols = []
        for imt in self.IMT_TYPES_ORDERED.keys():
            if (instrumented and 'MMI' in imt) or \
               (not instrumented and 'MMI' not in imt):
                continue
            new_cols.append(imt)
            ncols_amp += 1
        if ncols_amp > 0:
            app = np.empty((np.shape(station_rows)[0], ncols_amp))
            app[:] = np.nan

        col_dict = dict(zip(new_cols, range(ncols_amp)))
        id_dict = dict(zip(df['id'], range(nstation_rows)))

        #
        # Get all of the unflagged amps with the proper orientation
        #
        self.cursor.execute(
                'SELECT a.amp, a.imt_id, a.station_id FROM '
                'amp a, station s WHERE a.flag = "0" AND s.id = a.station_id '
                'AND s.instrumented = %d AND a.orientation NOT IN ("Z", "U") '
                'AND a.amp IS NOT NULL' % (instrumented)
            )
        amp_rows = self.cursor.fetchall()

        #
        # Go through all the amps and put them into the data frame
        #
        for this_row in amp_rows:
            imt_id = this_row[1]
            if imt_id not in self.IMT_LOOKUP:
                continue
            #
            # Set the cell to the peak amp
            #
            imt = self.IMT_LOOKUP[imt_id]
            rowidx = id_dict[this_row[2]]
            cval = app[rowidx, col_dict[imt]]
            amp = this_row[0]
            if np.isnan(cval) or (cval < amp):
                app[rowidx, col_dict[imt]] = amp

        if ncols_amp > 0:
            df2 = pd.DataFrame(app, columns=new_cols)
            df = pd.concat([df, df2], axis=1)

        del amp_rows, new_cols, col_dict, id_dict

        df['name'] = df.network.map(str) + '.' + df.code.map(str)

        if sort is True:
            if pd.__version__ >= '0.17.0':
                df = df.sort_values('name')
            else:
                df = df.sort('name')

        return df


    def getIMTset(self):
        """
        Return a set of IMTs in the StationList.

        Returns:
            Set of IMTs
        """

        return self.imtset.copy()


    @staticmethod
    def _getOrientation(orig_channel):
        """
        Return a character representing the orientation of a channel.

        Args:
            orig_channel (string):
                String representing the seed channel (e.g. 'HNZ'). The
                final character is assumed to be the (uppercase) orientation.

        Returns:
            Character representing the channel orientation. One of 'N',
            'E', 'Z', 'H' (for horizontal), or 'U' (for unknown).
        """
        if orig_channel[-1] in ('N', 'E', 'Z'):
            orientation = orig_channel[-1]
        elif orig_channel[-1] == "K":   # Channel is "UNK"; assume horizontal
            orientation = 'H'
        else:
            orientation = 'U'  # this is unknown
        return orientation


    @staticmethod
    def _getGroundMotions(comp, imt_translate):
        """
        Get a dictionary of peak ground motions (values and flags).
        Output keys are one of: [pga,pgv,psa03,psa10,psa30]
        Even if flags are not specified in the input, they will
        be guaranteed to at least have a flag of '0'.
        """
        pgmdict = {}
        imtset = set()
        for pgm in comp.childNodes:
            if pgm.nodeName == '#text':
                continue
            key = pgm.nodeName
            if key == 'acc':
                key = 'pga'
            elif key == 'vel':
                key = 'pgv'
            if key not in imt_translate:
                if 'pga' in key:
                    new_key = 'PGA'
                elif 'pgv' in key:
                    new_key = 'PGV'
                elif 'mmi' in key:
                    new_key = 'MMI'
                elif 'psa' in key:
                    pp = get_imt_period(key)
                    new_key = 'SA(' + str(pp) + ')'
                else:
                    raise ValueError('Unknown amp type in input: %s' % key)
#                    new_key = key
                imt_translate[key] = new_key
            else:
                new_key = imt_translate[key]
            key = new_key

            value = float(pgm.getAttribute('value'))
            if pgm.hasAttribute('flag'):
                flag = pgm.getAttribute('flag')
            else:
                flag = '0'
            pgmdict[key] = {'value': value, 'flag': flag}
            imtset.add(key)
        return pgmdict, imtset


    @staticmethod
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


    @staticmethod
    def _filter_station(xmlfile):
        """
        Filter individual xmlfile into a stationdict data structure.

        Args:
            xmlfile (string):
                Path to ShakeMap XML input file (or file-like object)
                containing station data.

        Returns:
            stationdict data structure
        """
        stationdict = OrderedDict()
        imt_translate = {}
        imtset = set()
        dom = minidom.parse(xmlfile)
        for root in dom.childNodes:
            if not isinstance(root, minidom.DocumentType):
                break
        stations = root.getElementsByTagName('station')
        for station in stations:
            code = station.getAttribute('code')
            attributes = StationList._getStationAttributes(station)
            comps = station.getElementsByTagName('comp')
            if code in stationdict:
                compdict = stationdict[code]
            else:
                compdict = {}
            for comp in comps:
                compname = comp.getAttribute('name')
                if 'Intensity Questionnaire' in str(compname):
                    compdict['mmi'] = {}
                    continue
                tpgmdict, ims = StationList._getGroundMotions(comp, imt_translate)
                if compname in compdict:
                    pgmdict = compdict[compname]
                else:
                    pgmdict = {}
                pgmdict.update(tpgmdict)
                # copy the VALUES, not REFERENCES, of the component list into
                # our growing dictionary
                compdict[compname] = copy.deepcopy(pgmdict)
                imtset |= ims
            if 'intensity' in attributes:
                compdict['mmi']['MMI'] = {
                    'value': attributes['intensity'], 'flag': '0'}
                imtset.add('MMI')
            stationdict[code] = (attributes, copy.deepcopy(compdict))
        dom.unlink()
        return stationdict, imtset


    @staticmethod
    def _createTables(db, cursor, imtset):
        """
        Build the database tables.
        """
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

        # IMT types
        rows = []
        for imt_id, imt_type in enumerate(imtset):
            rows.append((imt_id, imt_type))
        cursor.executemany('INSERT INTO imt (id, type) VALUES (?, ?)',
                rows)

        # Soil types
        for soiltype, sid in SOIL_TYPES.items():
            cursor.execute(
                    'INSERT INTO soiltype (soil_type, id) VALUES '
                    '("%s", %d)' % (soiltype, sid)
                )
        db.commit()

def get_imt_period(imt):

    p = re.search('(?<=psa)\d+', imt)
    return float(p.group(0)[:-1] + '.' + p.group(0)[-1])

def imt_sort(key):

    if isinstance(key, tuple):
        key = key[0]
    if key == 'MMI':
        return 0
    elif key == 'PGA':
        return 1
    elif key == 'PGV':
        return 2
    elif 'SA' in key:
        return 2 + float(key[3:-1])
    else:
        raise ValueError('Error in imt_sort, unknown key=%s' % str(key))
