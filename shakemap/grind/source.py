#!/usr/bin/env python

# stdlib imports
from xml.dom import minidom
import io
import sys
import copy

# third party
import pytz
import numpy as np
from openquake.hazardlib.geo import geodetic
from openquake.hazardlib.geo import point
from openquake.hazardlib.geo import Mesh
from openquake.hazardlib.geo.surface import multi
from openquake.hazardlib.geo.surface import planar
from openquake.hazardlib.geo import utils
from openquake.hazardlib.gsim import base, abrahamson_2014
from openquake.hazardlib.const import TRT
from shapely import geometry
from ..utils.timeutils import ShakeDateTime
from .ecef import latlon2ecef, ecef2latlon
from .fault import Fault

# local imports
from shakemap.utils.exception import ShakeMapException

REQUIRED_KEYS = ['lat', 'lon', 'depth', 'time', 'mag']
OPTIONAL_KEYS = ['id', 'timezone', 'locstring', 'type', 'created']
DEG2RAD = 180. / np.pi
RAKEDICT = {'SS': 0.0, 'NM': -90.0, 'RS': 90.0, 'ALL': None}
DEFAULT_STRIKE = 0.0
DEFAULT_DIP = 90.0
DEFAULT_WIDTH = 0.0
DEFAULT_ZTOR = 0.0


def read_event_file(eventxml):
    """
    Read event.xml file from disk, returning a dictionary of attributes.
    Input XML format looks like this:
    ******
    <earthquake id="2008ryan" lat="30.9858" lon="103.3639" mag="7.9" year="2008" month="05" day="12" hour="06" minute="28" second="01" timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA" created="1211173621" otime="1210573681" type="" />
    ******
    :param eventxml:
        Path to event XML file OR file-like object
    :returns:
       Dictionary with keys as indicated above in earthquake element attributes.
    """
    # fill in default values for mechanism, rake and dip
    # these may be overriden by values in event.xml, source.txt, or by values
    # passed in after reading input data.
    event = {'mech': 'ALL', 'rake': 45.0, 'dip': 90.0}
    if isinstance(eventxml, str):
        root = minidom.parse(eventxml)
    else:
        data = eventxml.read()
        root = minidom.parseString(data)

    eq = root.getElementsByTagName('earthquake')[0]
    event['lat'] = float(eq.getAttribute('lat'))
    event['lon'] = float(eq.getAttribute('lon'))
    event['depth'] = float(eq.getAttribute('depth'))
    event['mag'] = float(eq.getAttribute('mag'))
    year = int(eq.getAttribute('year'))
    month = int(eq.getAttribute('month'))
    day = int(eq.getAttribute('day'))
    hour = int(eq.getAttribute('hour'))
    minute = int(eq.getAttribute('minute'))
    second = float(eq.getAttribute('second'))
    msec = int((second - int(second)) * 1e6)
    second = int(second)

    # TODO: Handle:
    # 1) timezones other than UTC (use pytz)
    # 2) times < year 1900 (strftime fails) - subclass datetime, use Obspy utcdatetime
    # 3) daylight savings time (use pytz)

    # read in optional float fields
    for key in ['mech', 'rake', 'dip']:
        if eq.hasAttribute(key):
            event[key] = float(eq.getAttribute(key))

    event['time'] = ShakeDateTime(year, month, day, hour, minute, second, msec)
    # read in all of the fields that we know to be strings
    for key in ['id', 'locstring', 'type', 'timezone']:
        if eq.hasAttribute(key):
            event[key] = eq.getAttribute(key)
    if 'created' in event:
        event['created'] = ShakeDateTime.utcfromtimestamp(
            int(eq.getAttribute('created')))

    root.unlink()
    return event


def read_source(sourcefile):
    """
    Read source.txt file, which has lines like key=value.
    :param sourcefile:
        Path to source.txt file OR file-like object
    :returns:
        Dictionary containing key/value pairs from file.
    """
    isFile = False
    if isinstance(sourcefile, str):
        f = open(sourcefile, 'rt')
    else:
        isFile = True
        f = sourcefile
    params = {}
    for line in f.readlines():
        if line.startswith('#'):
            continue
        parts = line.split('=')
        key = parts[0].strip()
        value = parts[1].strip()
        # see if this is some sort of number
        try:
            value = int(value)
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                pass
        params[key] = value
    if not isFile:
        f.close()
    return params


class Source(object):
    """
    Encapsulate everything we need to know about an earthquake source.
    """

    def __init__(self, event, fault=None):
        """
        Construct a Source object.
        :param event:
            dictionary of values (see read_event_file() and read_source())
        :param fault:
            a Fault object
        :returns:
            Source object.
        """
        missing = []
        for req in REQUIRED_KEYS:
            if req not in list(event.keys()):
                missing.append(req)
        if len(missing):
            raise NameError('Input event dictionary is missing the following '
                            'required keys: "%s"' % (','.join(missing)))
        self._fault = fault
        self._event_dict = event.copy()

    @classmethod
    def readFromFile(cls, eventxmlfile, faultfile=None, sourcefile=None):
        """
        Class method to create a Source object by specifying an event.xml file,
        a fault file, and a source.txt file.
        :param eventxmlfile:
            Event xml file (see read_event_file())
        :param faultfile:
            Fault text file (see fault.py)
        :param sourcefile:
            source.txt file (see read_source())
        :returns:
            Source object
        """
        event = read_event_file(eventxmlfile)
        if faultfile is not None:
            fault = Fault.readFaultFile(faultfile)
        else:
            fault = None
        params = None
        if sourcefile is not None:
            params = read_source(sourcefile)
            event.update(params)
        return cls(event, fault=fault)

    def getRuptureContext(self, gmpelist):
        """
        Return a GEM RuptureContext suitable for any GMPE (except for those requiring hypo_loc).
        https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py
        :param gmpelist:
          Sequence of hazardlib GMPE objects.
        :raises:
          ShakeMapException when a GMPE requiring 'hypo_loc' is passed in.
        :returns:
            RuptureContext object with all known parameters filled in.
        If Source does not contain a Fault, strike, dip, ztor, and width will be
        filled with default values. Rake may not be known, or may be estimated
        from a focal mechanism.
        """
        # rupturecontext constructor inputs:
        # 'mag', 'strike', 'dip', 'rake', 'ztor', 'hypo_lon', 'hypo_lat',
        # 'hypo_depth', 'width', 'hypo_loc'
        for gmpe in gmpelist:
            reqset = gmpe.REQUIRES_RUPTURE_PARAMETERS
            if 'hypo_loc' in reqset:
                raise ShakeMapException(
                    'Rupture parameter "hypo_loc" is not supported!')

        rup = base.RuptureContext()
        rup.mag = self.getEventParam('mag')
        if self._fault is not None:
            rup.strike = self._fault.getStrike()
            rup.dip = self._fault.getDip()
            rup.ztor = self._fault.getTopOfRupture()
            rup.width = self._fault.getWidth()
        else:
            rup.strike = DEFAULT_STRIKE
            rup.dip = self.getEventParam('dip')
            rup.ztor = DEFAULT_ZTOR
            rup.width = DEFAULT_WIDTH

        if 'rake' in self._event_dict:
            rup.rake = self.getEventParam('rake')
        elif 'mech' in self._event_dict:
            mech = self._event_dict['mech']
            rup.rake = RAKEDICT[mech]
        else:
            rup.rake = RAKEDICT['ALL']

        rup.hypo_lat = self.getEventParam('lat')
        rup.hypo_lon = self.getEventParam('lon')
        rup.hypo_depth = self.getEventParam('depth')
        return rup

    def getFault(self):
        """
        Return fault instance.
        :returns:
           Fault instance
        """
        return copy.deepcopy(self._fault)

    def getHypo(self):
        """
        Return hypocenter point.
        :returns:
           Hypocenter point.
        """
        hyplat = self.getEventParam('lat')
        hyplon = self.getEventParam('lon')
        hypdepth = self.getEventParam('depth')
        hyppoint = point.Point(hyplon, hyplat, hypdepth)
        return hyppoint

    def setTectonicRegion(self, region):
        """
        Set tectonic region.
        :param region:
            TRT object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/const.py)
        """
        if not isinstance(region, TRT):
            raise ValueError(
                'Input tectonic region must be of type openquake.hazardlib.const.TRT')
        self._tectonic_region = copy.deepcopy(region)

    def getTectonicRegion(self):
        """
        Get tectonic region.
        :returns:
            TRT object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/const.py)
        """
        return copy.deepcopy(self.TectonicRegion)

    def getEventDict(self):
        """
        Return the event dictionary.
        :returns:
           Copy of dictionary of event keys/values
        """
        return self._event_dict.copy()

    def setEventParam(self, key, value):
        """
        Set a parameter in the internal event dictionary
        :param key:
            string key
        :param value:
            value (any object type)
        """
        self._event_dict[key] = value

    def setMechanism(self, mech, rake=None, dip=None):
        """
        Set the earthquake mechanism manually (overriding any values read in from event.xml or source.txt.
        If rake and dip are not specified, they will be assigned by mechanism as follows:
        Mech     Rake    Dip
        ----------------------
        RS         90     40
        NM        -90     50
        SS          0     90
        ALL        45     90

        :param mech:
          String - one of 'RS' (reverse), 'NM' (normal), 'SS' (strike slip), or 'ALL' (unknown).
        :param rake:
          Value between -360 and 360 degrees. If set, will override default value for mechanism (see table above).
        :param dip:
          Value betweeen 0 and 90 degrees. If set, will override default value for mechanism (see table above).
          (Value will be converted to range between -180 and 180 degrees.)
        """
        mechs = {'RS': {'rake': 90.0, 'dip': 40.0},
                 'NM': {'rake': -90.0, 'dip': 50.0},
                 'SS': {'rake': 0.0, 'dip': 90.0},
                 'ALL': {'rake': 45.0, 'dip': 90.0}}
        if mech not in list(mechs.keys()):
            raise ShakeMapException(
                'Mechanism must be one of: %s' % str(
                    list(
                        mechs.keys())))
        if dip is not None:
            if dip < 0 or dip > 90:
                raise ShakeMapException('Dip must be in range 0-90 degrees.')
            if dip < 0 or dip > 90:
                raise ShakeMapException('Dip must be in range 0-90 degrees.')
        else:
            dip = mechs[mech]['dip']

        if rake is not None:
            if rake < -180:
                rake += 360
            if rake > 180:
                rake -= 360
            if rake < -180 or rake > 180:
                raise ShakeMapException(
                    'Rake must be transformable to be in range -180 to 180 degrees.')
        else:
            rake = mechs[mech]['rake']

        self.setEventParam('dip', dip)
        self.setEventParam('rake', rake)
        self.setEventParam('mech', mech)

    def getEventParam(self, key):
        """
        Get a parameter from the internal event dictionary
        :param key:
            string key
        :returns:
            value (any object type)
        """
        if key not in list(self._event_dict.keys()):
            raise NameError('Key "%s" not found in event dictionary' % key)
        # rather than rais an error, can we return None?
        # e.g. self._event_dict.get(key)
        return copy.deepcopy(self._event_dict[key])

    def setFaultReference(self, reference):
        """
        Set the citeable reference for the fault
        :param reference:
            string citeable reference
        """
        self._fault.setReference(reference)

    def getFaultReference(self):
        """
        Get the citeable reference for the fault
        :returns:
            string citeable reference
        """
        return self._fault.getReference()

    def getQuadrilaterals(self):
        """
        Get the list of quadrilaterals defining the fault.
        :returns:
            List of quadrilateral tuples (4 Point objects each)
        """
        return self._fault.getQuadrilaterals()


def rake_to_mech(rake):
    mech = 'ALL'
    if (rake >= -180 and rake <= -150) or \
       (rake >= -30  and rake <= 30) or \
       (rake >= 150 and rake <= 180):
        mech = 'SS'
    if rake >= -120 and rake <= -60:
        mech = 'NM'
    if rake >= 60 and rake <= 120:
        mech = 'RS'
    return mech
