#!/usr/bin/env python

#stdlib imports
import sqlite3
from xml.dom import minidom
import os.path
import sys
import copy
import time

#third party imports
import pandas as pd
import numpy as np
from openquake.hazardlib import imt as GEM_IMT

#local imports
from shakemap.gmice.wgrw12 import WGRW12
from .distance import get_distance

TABLES = {'station':
          {'id':'integer primary key',
           'network':'str',
           'code':'str',
           'name':'str',
           'lat':'float',
           'lon':'float',
           'elev':'float',
           'repi':'float',
           'rhypo':'float',
           'rrup':'float',
           'rjb':'float',
           'vs30':'float',
           'instrumented':'int'},
          'imt':
           {'id':'integer primary key',
            'station_id':'int',
            'imt_type':'str'},
          'amp':
           {'id':'integer primary key',
            'station_id':'int',
            'imt_id':'int',
            'original_channel':'str',
            'orientation':'str',
            'amp':'float',
            'uncertainty':'float',
            'flag':'str'},
          'soiltype':
           {'id':'integer primary key',
            'soil_type':'str'},
          'predicted':
           {'id':'integer primary key',
            'soiltype_id':'int',
            'imt_id':'int',
            'station_id':'int',
            'amp':'float',
            'uncertainty':'float'}
            }

IMT_TYPES = ['mmi','pga','pgv','psa03','psa10','psa30',
             'pga_mmi','pgv_mmi','psa03_mmi','psa10_mmi','psa30_mmi',
             'mmi_pga','mmi_pgv','mmi_psa03','mmi_psa10','mmi_psa30']

#dictionary of our imt type strings to the kind that GEM needs to create IMT objects.
IMT_MAP = {'mmi':'MMI',
           'pga':'PGA',
           'pgv':'PGV',
           'psa03':'SA(0.3)',
           'psa10':'SA(1.0)',
           'psa30':'SA(3.0)'}

def _getOrientation(orig_channel):
    if orig_channel[-1] in ('N','E','Z'):
        orientation = orig_channel[-1]
    else:
        orientation = 'U' #this is unknown
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
        pgmdict[key] = {'value':value,'flag':flag}
    return pgmdict


def _getStationAttributes(station):
    """
    Get a dictionary of the station attributes
    """
    attrdict = {}
    for attr in list(station.attributes.items()):
        key = attr[0]
        value = attr[1]
        #is this value a float or str?
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
        if not isinstance(root,minidom.DocumentType):
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
            #copy the VALUES, not REFERENCES, of the component list into our growing dictionary
            compdict[compname] = copy.deepcopy(pgmdict)
        if 'intensity' in list(attributes.keys()):
            compdict[compname]['mmi'] = {'value':attributes['intensity'],'flag':'0'}
        stationdict[code] = (attributes,copy.deepcopy(compdict))
    dom.unlink()
    return stationdict

def _createTables(db,cursor):
    for table in TABLES.keys():
        sql = 'CREATE TABLE %s (' % table
        nuggets = []
        for column,ctype in TABLES[table].items():
            nuggets.append('%s %s' % (column,ctype))
        sql += ','.join(nuggets) + ')'
        cursor.execute(sql)
        db.commit()
    #IMT types are either observed (first row here)
    #derived MMI (second row)
    #or derived PGM (third row)
    
    for imt_type in IMT_TYPES:
        imtquery = 'INSERT INTO imt (imt_type) VALUES ("%s")' % imt_type
        cursor.execute(imtquery)
        db.commit()
    soiltypes = ['rock','soil']
    for soiltype in soiltypes:
        soilquery = 'INSERT INTO soiltype (soil_type) VALUES ("%s")' % soiltype
        cursor.execute(soilquery)
        db.commit()

class StationList(object):
    def __init__(self,dbfile):
        self.db = sqlite3.connect(dbfile)
        self.cursor = self.db.cursor()

    def __len__(self):
        squery = 'SELECT count(*) FROM station'
        self.cursor.execute(squery)
        return self.cursor.fetchone()[0]
        
    @classmethod
    def loadFromDict(cls,stationdictlist,dbfile):
        do_create = False
        if not os.path.isfile(dbfile):
            do_create = True
        dbfile = dbfile
        db = sqlite3.connect(dbfile)
        cursor = db.cursor()
        if do_create:
            #create the tables we want
            _createTables(db,cursor)

        for stationdict in stationdictlist:
            for key,station_tpl in stationdict.items():
                station_attributes,comp_dict = station_tpl
                lat = station_attributes['lat']
                lon = station_attributes['lon']
                network = station_attributes['netid']
                code = key
                if key.startswith(network):
                    code = key.replace(network+'.','')
                name = station_attributes['name']
                #elevation?
                ciim_tuple = ('dyfi','mmi','intensity','ciim')
                instrumented = int(station_attributes['netid'].lower() not in ciim_tuple) #????
                fmt = 'INSERT INTO station (network,code,name,lat,lon,instrumented) VALUES ("%s","%s","%s",%.4f,%.4f,%i)'
                station_query = fmt % (network,code,name,lat,lon,instrumented)
                cursor.execute(station_query)
                db.commit()
                sid = cursor.lastrowid
                for original_channel,pgm_dict in comp_dict.items():
                    orientation = _getOrientation(original_channel)
                    for imt_type,imt_dict in pgm_dict.items():
                        imtquery = 'SELECT id FROM imt WHERE imt_type == "%s"' % imt_type
                        cursor.execute(imtquery)
                        imtid = cursor.fetchone()[0]
                        amp = imt_dict['value']
                        if np.isnan(amp):
                            amp = 'NULL'
                            fmt = 'INSERT INTO amp (station_id,imt_id,original_channel,orientation,amp,flag) VALUES (%i,%i,"%s","%s","%s","%s")'
                        else:
                            fmt = 'INSERT INTO amp (station_id,imt_id,original_channel,orientation,amp,flag) VALUES (%i,%i,"%s","%s",%.4f,"%s")'
                        flag = imt_dict['flag']
                        
                        ampquery = fmt % (sid,imtid,original_channel,orientation,amp,flag)
                        try:
                            cursor.execute(ampquery)
                        except:
                            x = 1
                        db.commit()
            
        cursor.close()
        db.close()
        return cls(dbfile)

    def fillTables(self,source):
        """Populate tables with derived MMI/PGM values and distances.
        :param source:
          ShakeMap Source object.
        """
        gmice = WGRW12()
        #find all of the instrumented stations
        stationquery = 'SELECT id,lat,lon,code,network FROM station where instrumented = 1'
        self.cursor.execute(stationquery)
        rows = self.cursor.fetchall()
        emag = source.getEventDict()['mag']
        distances = ['rhypo','repi','rjb','rrup']

        #pre-fetch all of the IMT ids before looping
        imtdict = {}
        for imt in IMT_TYPES:
            imtquery = 'SELECT id FROM imt WHERE imt_type = "%s"' % imt
            self.cursor.execute(imtquery)
            try:
                imtdict[imt] = self.cursor.fetchone()[0]
            except:
                x = 1
            
        for row in rows:
            sid,lat,lon,code,network = row
            #calculate all distance types
            ddict = get_distance(distances,np.array([lat]),np.array([lon]),np.array([0]),source)
            values = []
            for d in distances:
                values.append(ddict[d][0])
            values.append(sid)
            values = tuple(values)
            station_update = 'UPDATE station set rhypo=%.2f,repi=%.2f,rjb=%.2f,rrup=%.2f WHERE id=%i' % values
            self.cursor.execute(station_update)
            self.db.commit()

            #calculate all derived mmi values
            for imt,imtid in imtdict.items():
                if imt.endswith('_mmi') or imt.startswith('mmi'):
                    continue
                #what distance measure to use here?
                ampquery = 'SELECT amp FROM amp WHERE station_id=%i AND imt_id=%i' % (sid,imtid)
                self.cursor.execute(ampquery)
                imtvalue = self.cursor.fetchone()[0]
                gemimt = GEM_IMT.from_string(IMT_MAP[imt])
                dmmi = gmice.getMIfromGM(imtvalue,gemimt,dists=ddict['repi'][0],mag=emag)
                derived_mmi = imt+'_mmi'
                derived_imtid = imtdict[derived_mmi]
                self.cursor.execute('INSERT INTO amp (imt_id,amp,station_id,flag) VALUES (%i,%.2f,%i,"0")' % (derived_imtid,dmmi,sid))
                self.db.commit()

            #calculate all derived pgm values
            mmiquery = 'SELECT amp FROM amp WHERE station_id=%i AND imt_id=%i' % (sid,imtdict['mmi'])
            self.cursor.execute(mmiquery)
            mmivalue = self.cursor.fetchone()[0]
            for imt,imtid in imtdict.items():
                if not imt.startswith('mmi_'):
                    continue
                #what distance measure to use here?
                tmp,derived_imt = imt.split('_')
                gemimt = GEM_IMT.from_string(IMT_MAP[derived_imt])
                dpgm = gmice.getGMfromMI(mmivalue,gemimt,dists=ddict['repi'],mag=emag)
                self.cursor.execute('INSERT INTO amp (imt_id,amp,station_id,flag) VALUES (%i,%.2f,%i,"0")' % (imtid,dpgm,sid))
                self.db.commit()
        
    def __del__(self):
        self.cursor.close()
        self.db.close()

    @classmethod
    def loadFromXML(cls,xmlfiles,dbfile):
        stationdictlist = []
        for xmlfile in xmlfiles:
            stationdict = _filter_station(xmlfile)
            stationdictlist.append(stationdict)
        return cls.loadFromDict(stationdictlist,dbfile)

    def getInstrumentedStations(self):
        stationquery = 'SELECT id,lat,lon,code,network FROM station where instrumented = 1'
        self.cursor.execute(stationquery)
        rows = self.cursor.fetchall()
        columns = ['id','lat','lon','code','network']
        try:
            df = pd.DataFrame(rows,columns=columns)
        except Exception as e:
            x = 1
        imts = ['pga','pgv','psa03','psa10','psa30','pga_mmi','pgv_mmi','psa03_mmi','psa10_mmi','psa30_mmi']
        for imt in imts:
            df[imt] = np.nan
            df[imt+'_unc'] = np.nan
        rowidx = 0
        for row in rows:
            sid = row[0]
            code = row[3]
            for imt in imts:
                if code.find('LPK') > -1 and imt == 'pgv':
                    c = 1
                imtquery = 'SELECT id FROM imt WHERE imt_type = "%s"' % imt
                self.cursor.execute(imtquery)
                imtid = self.cursor.fetchone()[0]
                ampquery = 'SELECT amp,uncertainty FROM amp WHERE amp.flag = "0" AND imt_id = %i AND station_id=%i' % (imtid,sid)
                self.cursor.execute(ampquery)
                row = self.cursor.fetchone()
                if row is not None:
                    if row[0] is None:
                        amp = np.nan
                    else:
                        amp = row[0]
                    if row[1] is None:
                        unc = np.nan
                    else:
                        unc = row[0]
                    df.ix[rowidx,imt] = amp
                    df.ix[rowidx,imt+'_unc'] = unc

            rowidx += 1

        df['name'] = df.network.map(str) + '.' + df.code.map(str)
        del df['network']
        del df['code']
        newcols = ['name','lat','lon',]
        for imt in imts:
            newcols.append(imt)
            newcols.append(imt+'_unc')
        if pd.__version__ >= '0.17.0':
            df = df[newcols].sort_values('name')
        else:
            df = df[newcols].sort('name')
        return df

    def getMMIStations(self):
        stationquery = 'SELECT id,lat,lon,code,network FROM station where instrumented = 0'
        self.cursor.execute(stationquery)
        rows = self.cursor.fetchall()
        columns = ['id','lat','lon','code','network']
        df = pd.DataFrame(rows,columns=columns)
        df['mmi'] = np.nan
        df['mmi_unc'] = np.nan
        rowidx = 0
        for row in rows:
            sid = row[0]
            imtquery = 'SELECT id FROM imt WHERE imt_type = "mmi"'
            self.cursor.execute(imtquery)
            imtid = self.cursor.fetchone()[0]
            ampquery = 'SELECT amp,uncertainty FROM amp WHERE amp.flag = "0" AND imt_id = %i AND station_id=%i' % (imtid,sid)
            self.cursor.execute(ampquery)
            row = self.cursor.fetchone()
            if row is not None:
                df.ix[rowidx,'mmi'] = row[0]
                df.ix[rowidx,'mmi_unc'] = row[1]

            rowidx += 1

        df['name'] = df.network.map(str) + '.' + df.code.map(str)
        del df['network']
        del df['code']
        newcols = ['name','lat','lon','mmi','mmi_unc']
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
    stations = StationList.loadFromXML(dbfile,xmlfiles)
    t2 = time.time()
    print('%i stations loaded in %.2f seconds' % (len(stations),t2-t1))
    mmidf = stations.getMMIStations()
    t3 = time.time()
    print('%i intensity observations retrieved in %.2f seconds' % (len(mmidf),t3-t2))
    imtdf = stations.getInstrumentedStations()
    t4 = time.time()
    print('%i instrumental measurements retrieved in %.2f seconds' % (len(imtdf),t4-t3))
    
