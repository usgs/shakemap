#!/usr/bin/env python

#stdlib modules
import sys
import os.path
from xml.dom import minidom
import copy
import pprint
import StringIO

#third party modules
import numpy as np

#local modules
from neicio.tag import Tag

def _getStationAttributes(station):
    """
    Get a dictionary of the station attributes
    """
    attrdict = {}
    for attr in station.attributes.items():
        key = attr[0]
        value = attr[1]
        #is this value a float or str?
        try:
            value = float(value)
        except:
            pass
        attrdict[key] = value
    return attrdict

def getPeakValues(stationdict):
    """Return the peak ground motion values (horizontal channels), intensity, station code, and coordinates of station data.
    Inputs:
    
     * stationdict Dictionary returned by filterStations()
     
    Output:
    
     * Dictionary containing numpy arrays for keys: lat,lon,pga,pgv,psa03,psa10,psa30,mmi,code

     NB: If any of the peak ground motion fields are not defined, numpy.NaN is assigned for that value,
    ensuring each array is of the same length.
    """
    peakdict = {'lat':[],
                'lon':[],
                'pga':[],
                'pgv':[],
                'psa03':[],
                'psa10':[],
                'psa30':[],
                'mmi':[],
                'code':[]
                }
    peakkeys = ['pga','pgv','psa03','psa10','psa30']
    for stationcode,stationtuple in stationdict.iteritems():
        attributes,compdict = stationtuple
        #we need to check the components to see if the two horizontal channels are distinguishable from the vertical channel
        channelsum = 0
        channels = compdict.keys()
        for channel in compdict.keys():
            if channel.endswith('1') or channel.endswith('2') or channel.endswith('3'):
                channelsum += 1
        if channelsum == 3:
            continue
        if attributes.has_key('intensity'):
            peakdict['mmi'].append(attributes['intensity'])
        else:
            peakdict['mmi'].append(np.nan)
        peakdict['lat'].append(attributes['lat'])
        peakdict['lon'].append(attributes['lat'])
        peakdict['code'].append(stationcode)
        #make sure we get the peak value from 
        for key in peakkeys:
            value = -1
            for channel,pgmdict in compdict.iteritems():
                if channel.lower().endswith('z'): #is vertical channel guaranteed to end with 'Z'????
                    continue
                if pgmdict.has_key(key) and pgmdict[key]['value'] > value and pgmdict[key]['flag'] == '0':
                    value = pgmdict[key]['value']
            if value == -1:
                value = np.nan
            peakdict[key].append(value)
    for key in peakdict.keys():
        peakdict[key] = np.array(peakdict[key])
    return peakdict
            
def writeStations(stationdict):
    """Write station dictionary structure to an XML string.
    Inputs:

     * stationdict Dictionary structure as returned by filterStations

    Outputs:
     * XML string suitable for writing to a file
    """
    earthquakeTag = Tag('earthquake-data') #top level
    for stationcode,stationtuple in stationdict.iteritems():
        attributes,compdict = stationtuple
        stationTag = Tag('station',attributes=attributes)
        for compname,pgmdict in compdict.iteritems():
            compTag = Tag('comp',attributes={'name':compname})
            for pgm,pgmvaldict in pgmdict.iteritems():
                pgmTag = Tag(pgm,attributes=pgmvaldict)
                compTag.addChild(pgmTag)
            stationTag.addChild(compTag)
        earthquakeTag.addChild(stationTag)
    xmlstr = earthquakeTag.renderToXML(ntabs=0)
    return xmlstr

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
        if pgm.hasAttribute('flag'):
            flag = pgm.getAttribute('flag')
        else:
            flag = '0'
        pgmdict[key] = {'value':value,'flag':flag}
    return pgmdict

def filterStations(xmlfiles):
    """Read in N station data files in XML format and condense them into one data structure.
    Input:
    
     * xmlfiles - Sequence of xml files containing station data
     
    Output:
     * stationdict, a dictionary (keys are station codes) of unique station two-tuples.  
       First element contains a dictionary of station attributes. The second element contains:
    
       * a dictionary of components (keys are component names), values are:
       
         * a dictionary of peak ground motions where keys are pga, pgv, etc. and values are a dictionary with keys {value,flag}.

    Example:
    Input XML file that looks like this:
    ####################################################################################
    <shakemap-data code_version="3.5" map_version="4">
      <earthquake id="us10001rvu" lat="-4.7599" lon="152.5556" mag="7.5" year="2015" month="3" day="29" hour="23" minute="48" second="31" timezone="GMT" depth="40" network="us" locstring="NEW BRITAIN REGION, PAPUA NEW GUINEA" created="1427728491" />
        <stationlist created="1427728491">
          <station code="IU.PMG" name="Port Moresby, New Guinea" insttype="UNK" lat="-9.4047" lon="147.1597" dist="731.845032" source="NEIC" netid="IU" commtype="UNK" loc="" intensity="1.5">
            <comp name="HNZ">
              <pga value="0.0207" flag="0" />
              <pgv value="0.0397" flag="0" />
              <psa03 value="0.0420" flag="0" />
              <psa10 value="0.0327" flag="0" />
              <psa30 value="0.0132" flag="0" />
            </comp>
            <comp name="HN2">
              <pga value="0.0218" flag="0" />
              <pgv value="0.0565" flag="0" />
              <psa03 value="0.0321" flag="0" />
              <psa10 value="0.0472" flag="0" />
              <psa30 value="0.0179" flag="0" />
            </comp>
            <comp name="HN1">
              <pga value="0.0227" flag="0" />
              <pgv value="0.0385" flag="0" />
              <psa03 value="0.0399" flag="0" />
              <psa10 value="0.0506" flag="0" />
              <psa30 value="0.0162" flag="0" />
            </comp>
          </station>
        </stationlist>
      </shakemap-data>
    ####################################################################################
    will be parsed into a structure that looks like this:
    cmpdict = { u'IU.PMG': ( { u'code': u'IU.PMG',
                 u'commtype': u'UNK',
                 u'dist': 731.845032,
                 u'insttype': u'UNK',
                 u'intensity': 1.5,
                 u'lat': -9.4047,
                 u'loc': u'',
                 u'lon': 147.1597,
                 u'name': u'Port Moresby, New Guinea',
                 u'netid': u'IU',
                 u'source': u'NEIC'},
               { u'HN1': { u'pga': { 'flag': True, 'value': 0.0227},
                           u'pgv': { 'flag': True, 'value': 0.0385},
                           u'psa03': { 'flag': True, 'value': 0.0399},
                           u'psa10': { 'flag': True, 'value': 0.0506},
                           u'psa30': { 'flag': True, 'value': 0.0162}},
                 u'HN2': { u'pga': { 'flag': True, 'value': 0.0218},
                           u'pgv': { 'flag': True, 'value': 0.0565},
                           u'psa03': { 'flag': True, 'value': 0.0321},
                           u'psa10': { 'flag': True, 'value': 0.0472},
                           u'psa30': { 'flag': True, 'value': 0.0179}},
                 u'HNZ': { u'pga': { 'flag': True, 'value': 0.0207},
                           u'pgv': { 'flag': True, 'value': 0.0397},
                           u'psa03': { 'flag': True, 'value': 0.042},
                           u'psa10': { 'flag': True, 'value': 0.0327},
                           u'psa30': { 'flag': True, 'value': 0.0132}}})}
    """
    stationdict = {}
    for xmlfile in xmlfiles:
        if not os.path.isfile(xmlfile):
            raise Exception('Input xml file "%s" does not exist.' % xmlfile)
        tstationdict = filterStation(xmlfile)
        stationdict.update(tstationdict)        
    return stationdict

def filterStation(xmlfile):
    """
    Filter individual xmlfile into a stationdict data structure.
    Inputs:

     * xmlfile xml file (or file-like object) containing station data

    Outputs:

     * stationdict Data structure as returned by filterStations()
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
        if code in stationdict.keys():
            compdict = stationdict[code]
        else:
            compdict = {}
        for comp in comps:
            compname = comp.getAttribute('name')
            tpgmdict = _getGroundMotions(comp)
            if compname in compdict.keys():
                pgmdict = compdict[compname]
            else:
                pgmdict = {}
            pgmdict.update(tpgmdict)
            #copy the VALUES, not REFERENCES, of the component list into our growing dictionary
            compdict[compname] = copy.deepcopy(pgmdict)
        stationdict[code] = (attributes,copy.deepcopy(compdict))
    dom.unlink()
    return stationdict

def _main(xmlfiles):
    stationdict = filterStations(txmlfiles)
    peakdict = getPeakValues(stationdict)
    printer = pprint.PrettyPrinter(indent=2)
    printer.pprint(stationdict)
    xmlstr = writeStations(stationdict)
    print xmlstr

def _test():
    xmlstr = """<shakemap-data code_version="3.5" map_version="4">
      <earthquake id="us10001rvu" lat="-4.7599" lon="152.5556" mag="7.5" year="2015" month="3" day="29" hour="23" minute="48" second="31" timezone="GMT" depth="40" network="us" locstring="NEW BRITAIN REGION, PAPUA NEW GUINEA" created="1427728491" />
        <stationlist created="1427728491">
          <station code="IU.PMG" name="Port Moresby, New Guinea" insttype="UNK" lat="-9.4047" lon="147.1597" dist="731.845032" source="NEIC" netid="IU" commtype="UNK" loc="" intensity="1.5">
            <comp name="HNZ">
              <pga value="0.0207" flag="0" />
              <pgv value="0.0397" flag="0" />
              <psa03 value="0.0420" flag="0" />
              <psa10 value="0.0327" flag="0" />
              <psa30 value="0.0132" flag="0" />
            </comp>
            <comp name="HN2">
              <pga value="0.0218" flag="0" />
              <pgv value="0.0565" flag="0" />
              <psa03 value="0.0321" flag="0" />
              <psa10 value="0.0472" flag="0" />
              <psa30 value="0.0179" flag="0" />
            </comp>
            <comp name="HN1">
              <pga value="0.0227" flag="0" />
              <pgv value="0.0385" flag="0" />
              <psa03 value="0.0399" flag="0" />
              <psa10 value="0.0506" flag="0" />
              <psa30 value="0.0162" flag="0" />
            </comp>
          </station>
        </stationlist>
      </shakemap-data>"""
    cmpdict = { u'IU.PMG': ( { u'code': u'IU.PMG',
                 u'commtype': u'UNK',
                 u'dist': 731.845032,
                 u'insttype': u'UNK',
                 u'intensity': 1.5,
                 u'lat': -9.4047,
                 u'loc': u'',
                 u'lon': 147.1597,
                 u'name': u'Port Moresby, New Guinea',
                 u'netid': u'IU',
                 u'source': u'NEIC'},
               { u'HN1': { u'pga': { 'flag': True, 'value': 0.0227},
                           u'pgv': { 'flag': True, 'value': 0.0385},
                           u'psa03': { 'flag': True, 'value': 0.0399},
                           u'psa10': { 'flag': True, 'value': 0.0506},
                           u'psa30': { 'flag': True, 'value': 0.0162}},
                 u'HN2': { u'pga': { 'flag': True, 'value': 0.0218},
                           u'pgv': { 'flag': True, 'value': 0.0565},
                           u'psa03': { 'flag': True, 'value': 0.0321},
                           u'psa10': { 'flag': True, 'value': 0.0472},
                           u'psa30': { 'flag': True, 'value': 0.0179}},
                 u'HNZ': { u'pga': { 'flag': True, 'value': 0.0207},
                           u'pgv': { 'flag': True, 'value': 0.0397},
                           u'psa03': { 'flag': True, 'value': 0.042},
                           u'psa10': { 'flag': True, 'value': 0.0327},
                           u'psa30': { 'flag': True, 'value': 0.0132}}})}
    xmlfile = StringIO.StringIO(xmlstr)
    stationdict = filterStation(xmlfile)
    assert _compareDictionaries(stationdict,cmpdict)
    #printer = pprint.PrettyPrinter(indent=2)
    #printer.pprint(stationdict)

def _compareDictionaries(dict1,dict2):
    """Compare two stationdict dictionaries at each level testing for equality.
    """
    #TODO: MAKE THIS MORE COMPLETE!!!!
    d1keys = set(dict1.keys())
    d2keys = set(dict2.keys())
    if len(d1keys.difference(d2keys)):
        return False
    return True
    
    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        _test()
        sys.exit(0)
    else:
        txmlfiles = sys.argv[1:]
        _main(txmlfiles)
    
    
