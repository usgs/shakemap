# stdlib imports
import sqlite3
from xml.dom import minidom
import os.path
import copy
from collections import OrderedDict

# third party imports
import pandas as pd
import numpy as np
from openquake.hazardlib import imt as GEM_IMT
from openquake.hazardlib.gsim import base

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

#
# This is an ordered version of IMT_TYPES just to help make the tables
# a little more human readable.
#
IMT_TYPES_ORDERED = OrderedDict(sorted(IMT_TYPES.items(), key=lambda t: t[1]))

#
# This flips the key/value pairs in the IMT_TYPES dictionary
#
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

    #
    # These are the netid's that indicate MMI data
    #
    __CIIM_TUPLE = ('dyfi', 'mmi', 'intensity', 'ciim')

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
        self.db = sqlite3.connect(dbfile)
        self.cursor = self.db.cursor()


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
        self.cursor.close()
        self.db.close()


    @classmethod
    def fromXML(cls, xmlfiles, dbfile, source, sites, gmpe, ipe, gmice):
        """
        Create a StationList object by reading one or more ShakeMap XML 
        input files and populate database tables with derived MMI/PGM 
        values and distances. This is a convenience method that calls 
        :meth:`loadFromXML` and :meth:`fillTables`.

        Args:
            xmlfiles (sequence of strings):
                Sequence of ShakeMap XML input files to read.
            dbfile (string): 
                Path to a file into which to write the SQLite database.
            source:
                ShakeMap Source object containing information about the 
                origin and source of the earthquake.
            sites:
                ShakeMap Sites object containing grid of Vs30 values for the 
                region in question.
            gmpe:
                A GMPE object to use for predicting 
                ground motions at the station locations.
            ipe:
                A GMPE object to use for predicting 
                macroseismic intensities at the station locations.
            gmice:
                A GMICE object to use for converting ground motions to
                macroseismic intensities.

        Returns:
            :class:`StationList` object
        """
    
        this = cls.loadFromXML(xmlfiles, dbfile)
        this.fillTables(source, sites, gmpe, ipe, gmice)
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
        for xmlfile in xmlfiles:
            stationdict = cls._filter_station(xmlfile)
            stationdictlist.append(stationdict)
        return cls._loadFromDict(stationdictlist, dbfile)


    @classmethod
    def _loadFromDict(cls, stationdictlist, dbfile):
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
        do_create = False
        if not os.path.isfile(dbfile):
            do_create = True
        dbfile = dbfile
        db = sqlite3.connect(dbfile)
        cursor = db.cursor()
        if do_create:
            # create the tables we want
            cls._createTables(db, cursor)

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
                                   not in cls.__CIIM_TUPLE)
                station_rows.append((sid, network, code, name, lat, lon, 
                                     instrumented))

                for original_channel, pgm_dict in comp_dict.items():
                    orientation = cls._getOrientation(original_channel)
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


    def fillTables(self, source, sites, gmpe, ipe, gmice):
        """
        Populate database tables with derived MMI/PGM values and distances.
        This method should be called after :meth:`loadFromXML`.

        Args:
            source:
                ShakeMap Source object containing information about the 
                origin and source of the earthquake.
            sites:
                ShakeMap Sites object containing grid of Vs30 values for the 
                region in question.
            gmpe:
                A GMPE object to use for predicting 
                ground motions at the station locations.
            ipe:
                A GMPE object to use for predicting 
                macroseismic intensities at the station locations.
            gmice:
                A GMICE object to use for converting ground motions to
                macroseismic intensities.

        Returns:
            nothing
        """
        emag = source.getEventDict()['mag']
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
        ddict = get_distance(DISTANCES, lats, lons, depths, source)
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
        # Use the GMICE to get the estimated MMIs from the various 
        # ground motion parameters from the instrumented stations.
        # These values will fill the rows like 'mmi_from_pgm'
        #
        amp_rows = []
        for imt, imtid in IMT_TYPES_ORDERED.items():
            if imt.endswith('_mmi') or imt.startswith('mmi'):
                continue
            self.cursor.execute(
                    'SELECT a.original_channel, a.orientation, a.uncertainty, '
                    'a.amp, s.repi, s.id FROM amp a, station s '
                    'WHERE a.imt_id=%i AND a.station_id=s.id AND '
                    's.instrumented=1 AND a.flag=0' % (imtid)
                )
            rows = self.cursor.fetchall()

            nrows = len(rows)
            amps = np.empty((nrows))
            dists = np.empty((nrows))
            emags = np.zeros((nrows)) + emag
            for irow, row in enumerate(rows):
                amps[irow] = row[3]
                dists[irow] = row[4]

            gemimt = GEM_IMT.from_string(GEM_IMT_MAP[imt])
            dmmi = gmice.getMIfromGM(amps, gemimt, dists=dists, mag=emags)

            derived_imtid = IMT_TYPES["mmi_from_" + imt]

            for irow, row in enumerate(rows):
                amp_rows.append((row[5], derived_imtid, row[0],
                                 row[1], dmmi[irow], 0.0))

        self.cursor.executemany(
            'INSERT INTO amp (station_id, imt_id, original_channel, '
            'orientation, amp, uncertainty, flag) VALUES '
            '(?, ?, ?, ?, ?, ?, "0")', amp_rows)

        # 
        # Use the GMICE to get the estimated ground motions from the 
        # mmi from the non-instrumented stations.
        # These will appeara in columns named like 'pgm_from_mmi'
        #
        imtid = IMT_TYPES['mmi']
        self.cursor.execute(
                'SELECT a.uncertainty, a.amp, s.repi, s.id FROM amp a, '
                'station s WHERE a.imt_id=%i AND a.station_id=s.id AND '
                's.instrumented=0 AND a.flag=0' % (imtid)
            )
        rows = self.cursor.fetchall()

        nrows = len(rows)
        mmi = np.empty((nrows))
        dists = np.empty((nrows))
        emags = np.zeros((nrows)) + emag
        for irow, row in enumerate(rows):
            mmi[irow] = row[1]
            dists[irow] = row[2]

        amp_rows = []
        for imt, imtid in IMT_TYPES_ORDERED.items():
            if imt.endswith('_mmi') or imt.startswith('mmi'):
                continue

            gemimt = GEM_IMT.from_string(GEM_IMT_MAP[imt])
            dmmi = gmice.getGMfromMI(mmi, gemimt, dists=dists, mag=emags)

            derived_imtid = IMT_TYPES[imt + "_from_mmi"]

            for irow, row in enumerate(rows):
                amp_rows.append((row[3], derived_imtid, dmmi[irow], 0.0))

        self.cursor.executemany(
            'INSERT INTO amp (station_id, imt_id, amp, uncertainty, flag) '
            'VALUES (?, ?, ?, ?, "0")', amp_rows)
        self.db.commit()

        #
        # Use the GMPE to predict ground motions for each of the
        # IMTs at each station
        # Compute predictions on both rock and soil, and find the
        # site amplification factors for each IMT
        # Also store the Vs30 for each station
        #
        rx = source.getRuptureContext([gmpe])
        dx = base.DistancesContext()
        for method in DISTANCES:
            (dx.__dict__)[method] = ddict[method]
        lldict = { 'lats': lats, 'lons': lons }
        sx_soil = sites.getSitesContext(lldict)
        sx_rock = sites.getSitesContext(lldict, rock_vs30=760.0)

        pred_rows = []
        siteamp_rows = []
        vs30_rows = []
        for imt in BASE_IMTS:
            gemimt = GEM_IMT.from_string(GEM_IMT_MAP[imt])
            if imt == 'mmi':
                stddev_types = ipe.DEFINED_FOR_STANDARD_DEVIATION_TYPES
                pred_soil, pred_stdev = ipe.get_mean_and_stddevs(
                        sx_soil, rx, dx, gemimt, stddev_types)
            else:
                stddev_types = gmpe.DEFINED_FOR_STANDARD_DEVIATION_TYPES
                pred_soil, pred_stdev = gmpe.get_mean_and_stddevs(
                        sx_soil, rx, dx, gemimt, stddev_types)
            
            if imt == 'mmi':
                pred_rock, junk = ipe.get_mean_and_stddevs(
                        sx_rock, rx, dx, gemimt, stddev_types)
            else:
                pred_rock, junk = gmpe.get_mean_and_stddevs(
                        sx_rock, rx, dx, gemimt, stddev_types)
            
            amp_facts = pred_soil - pred_rock

            if len(stddev_types) == 3:
                for irow, row in enumerate(station_rows):
                    pred_rows.append(
                            (row[0], IMT_TYPES[imt], 
                            SOIL_TYPES['rock'], pred_rock[irow], 
                            pred_stdev[0][irow], pred_stdev[1][irow], 
                            pred_stdev[2][irow])
                        )
                    siteamp_rows.append(
                            (row[0], IMT_TYPES[imt], amp_facts[irow])
                        )
            else:
                for irow, row in enumerate(station_rows):
                    pred_rows.append(
                            (row[0], IMT_TYPES[imt], 
                            SOIL_TYPES['rock'], pred_rock[irow], 
                            pred_stdev[0][irow], 'NULL', 'NULL')
                        )
                    siteamp_rows.append(
                            (row[0], IMT_TYPES[imt], amp_facts[irow])
                        )

        for irow, row in enumerate(station_rows):
            vs30_rows.append((sx_soil.vs30[irow], row[0]))
            
        self.cursor.executemany(
            'INSERT INTO predicted (station_id, imt_id, soiltype_id, amp, '
            'uncertainty_total, uncertainty_inter, uncertainty_intra) '
            'VALUES (?, ?, ?, ?, ?, ?, ?)', pred_rows)
        self.cursor.executemany(
            'INSERT INTO siteamp (station_id, imt_id, amp_factor) '
            'VALUES (?, ?, ?)', siteamp_rows)
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
        'rrup', 'rx', 'ry', 'ry0', 'U', 'T', 'mmi', 'mmi_unc', 'pga', 'pga_unc',
        'pgv', 'pgv_unc', 'psa03', 'psa03_unc', 'psa10', 'psa10_unc', 'psa30',
        'psa30_unc', 'mmi_from_pga', 'mmi_from_pga_unc', 'mmi_from_pgv',
        'mmi_from_pgv_unc', 'mmi_from_psa03', 'mmi_from_psa03_unc',
        'mmi_from_psa10', 'mmi_from_psa10_unc', 'mmi_from_psa30',
        'mmi_from_psa30_unc', 'mmi_predicted', 'mmi_predicted_unc_total',
        'mmi_predicted_unc_intra', 'mmi_site_factor', 'pga_predicted',
        'pga_predicted_unc_total', 'pga_predicted_unc_intra', 'pga_site_factor',
        'pgv_predicted', 'pgv_predicted_unc_total', 'pgv_predicted_unc_intra',
        'pgv_site_factor', 'psa03_predicted', 'psa03_predicted_unc_total',
        'psa03_predicted_unc_intra', 'psa03_site_factor', 'psa10_predicted',
        'psa10_predicted_unc_total', 'psa10_predicted_unc_intra',
        'psa10_site_factor', 'psa30_predicted', 'psa30_predicted_unc_total',
        'psa30_predicted_unc_intra', 'psa30_site_factor', 'name'

        For the non-instrumented dataframe, the columns would be:

        'id', 'lat', 'lon', 'code', 'network', 'vs30', 'repi', 'rhypo', 'rjb',
        'rrup', 'rx', 'ry', 'ry0', 'U', 'T', 'mmi', 'mmi_unc', 'pga', 'pga_unc',
        'pgv', 'pgv_unc', 'psa03', 'psa03_unc', 'psa10', 'psa10_unc', 'psa30',
        'psa30_unc', 'pga_from_mmi', 'pga_from_mmi_unc', 'pgv_from_mmi',
        'pgv_from_mmi_unc', 'psa03_from_mmi', 'psa03_from_mmi_unc',
        'psa10_from_mmi', 'psa10_from_mmi_unc', 'psa30_from_mmi',
        'psa30_from_mmi_unc', 'mmi_predicted', 'mmi_predicted_unc_total',
        'mmi_predicted_unc_intra', 'mmi_site_factor', 'pga_predicted',
        'pga_predicted_unc_total', 'pga_predicted_unc_intra', 'pga_site_factor',
        'pgv_predicted', 'pgv_predicted_unc_total', 'pgv_predicted_unc_intra',
        'pgv_site_factor', 'psa03_predicted', 'psa03_predicted_unc_total',
        'psa03_predicted_unc_intra', 'psa03_site_factor', 'psa10_predicted',
        'psa10_predicted_unc_total', 'psa10_predicted_unc_intra',
        'psa10_site_factor', 'psa30_predicted', 'psa30_predicted_unc_total',
        'psa30_predicted_unc_intra', 'psa30_site_factor', 'name'

        The **name** column is **network** and **code** concatenated with a 
        period (".") between them.
        The **unc** in column names means uncertainty. All ground motion, 
        site, and uncertainty 
        units are natural log units. Distances are in km. The intra-event 
        uncertainty (**<param>_predicted_unc_intra**) will be **np.nan** if 
        the GMPE or IPE doesn't define 
        separate inter- and intra-event terms.
        

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

        ncols = 0
        new_cols = []
        for imt in IMT_TYPES_ORDERED.keys():
            if (instrumented and imt.endswith('_from_mmi')) or \
               (not instrumented and imt.startswith('mmi_from_')):
                continue
            new_cols.append(imt)
            new_cols.append(imt + '_unc')
            ncols += 2
        app = np.empty((np.shape(station_rows)[0], ncols))
        app[:] = np.nan

        col_dict = dict(zip(new_cols, range(ncols)))
        id_dict = dict(zip(df['id'], range(nstation_rows)))

        #
        # Get all of the unflagged instrumented amps with the proper
        # orientation
        #
        self.cursor.execute(
            'SELECT a.amp, a.uncertainty, a.imt_id, a.station_id FROM '\
            'amp a, station s WHERE a.flag = "0" AND s.id = a.station_id '\
            'AND s.instrumented = %d AND a.orientation NOT IN ("Z", "U") '\
            'AND a.amp IS NOT NULL' % (instrumented)
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

        del amp_rows, new_cols, col_dict, id_dict
        #
        # Get all of the predicted amps along with their uncertainties
        #
        nstation_rows = len(station_rows)
        ncols = 0
        pred_cols = []
        for imt in BASE_IMTS:
            pred_cols.append(imt + '_predicted')
            pred_cols.append(imt + '_predicted_unc_total')
            pred_cols.append(imt + '_predicted_unc_intra')
            pred_cols.append(imt + '_site_factor')
            ncols += 4
        pred_arr = np.empty((np.shape(station_rows)[0], ncols))
        pred_arr[:] = np.nan

        col_dict = dict(zip(pred_cols, range(ncols)))
        id_dict = dict(zip(df['id'], range(nstation_rows)))

        self.cursor.execute(
            'SELECT p.amp, p.uncertainty_total, p.uncertainty_intra, '\
                   'p.imt_id, p.station_id FROM predicted p, station s '\
                   'WHERE s.id = p.station_id AND s.instrumented = %d AND '\
                   'p.amp IS NOT NULL' % (instrumented)
            )
        pred_rows = self.cursor.fetchall()

        #
        # Go through all the predicted amps and put them into the data 
        # frame
        #
        for this_row in pred_rows:
            imt_id = this_row[3]
            if imt_id not in IMT_LOOKUP:
                continue
            imt = IMT_LOOKUP[imt_id]
            rowidx = id_dict[this_row[4]]
            pred_arr[rowidx, col_dict[imt + '_predicted']] = this_row[0]
            pred_arr[rowidx, col_dict[imt + '_predicted_unc_total']] = this_row[1]
            if this_row[2] != 'NULL':
                pred_arr[rowidx, col_dict[imt + '_predicted_unc_intra']] = this_row[2]

        del pred_rows

        #
        # Go through all the site amplicication factors and put them 
        # into the data frame
        #
        self.cursor.execute(
                'SELECT p.amp_factor, p.imt_id, p.station_id FROM siteamp p, '
                'station s WHERE s.id = p.station_id AND '
                's.instrumented = %d AND p.amp_factor IS NOT NULL' % 
                (instrumented)
            )
        siteamp_rows = self.cursor.fetchall()

        for this_row in siteamp_rows:
            imt_id = this_row[1]
            if imt_id not in IMT_LOOKUP:
                continue
            imt = IMT_LOOKUP[imt_id]
            rowidx = id_dict[this_row[2]]
            pred_arr[rowidx, col_dict[imt + '_site_factor']] = this_row[0]

        df3 = pd.DataFrame(pred_arr, columns=pred_cols)
        df = pd.concat([df, df3], axis=1)

        df['name'] = df.network.map(str) + '.' + df.code.map(str)

        if sort is True:
            if pd.__version__ >= '0.17.0':
                df = df.sort_values('name')
            else:
                df = df.sort('name')

        return df


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
            'E', 'Z', or 'U' (for unknown).
        """
        if orig_channel[-1] in ('N', 'E', 'Z'):
            orientation = orig_channel[-1]
        else:
            orientation = 'U'  # this is unknown
        return orientation


    @staticmethod
    def _getGroundMotions(comp):
        """
        Get a dictionary of peak ground motions (values and flags).  
        Output keys are one of: [pga,pgv,psa03,psa10,psa30]
        Even if flags are not specified in the input, they will 
        be guaranteed to at least have a flag of '0'.
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
            if pgm.hasAttribute('flag'):
                flag = pgm.getAttribute('flag')
            else:
                flag = '0'
            pgmdict[key] = {'value': value, 'flag': flag}
        return pgmdict


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
        stationdict = {}
        dom = minidom.parse(xmlfile)
        for root in dom.childNodes:
            if not isinstance(root, minidom.DocumentType):
                break
        stations = root.getElementsByTagName('station')
        for station in stations:
            code = station.getAttribute('code')
            attributes = StationList._getStationAttributes(station)
            comps = station.getElementsByTagName('comp')
            if code in list(stationdict.keys()):
                compdict = stationdict[code]
            else:
                compdict = {}
            for comp in comps:
                compname = comp.getAttribute('name')
                tpgmdict = StationList._getGroundMotions(comp)
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


    @staticmethod
    def _createTables(db, cursor):
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
