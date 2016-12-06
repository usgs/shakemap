#!/usr/bin/env python

# stdlib imports
from xml.dom import minidom
import copy
import time as time

# third party
import numpy as np
from openquake.hazardlib.geo import geodetic
from openquake.hazardlib.geo import point
from openquake.hazardlib.gsim import base
from openquake.hazardlib.const import TRT

from shakemap.utils.timeutils import ShakeDateTime
from shakemap.grind.rupture import read_rupture_file
from shakemap.utils.exception import ShakeMapException


REQUIRED_KEYS = ['lat', 'lon', 'depth', 'mag']

RAKEDICT = {'SS': 0.0, 'NM': -90.0, 'RS': 90.0, 'ALL': None}
DEFAULT_MECH = 'ALL'
DEFAULT_STRIKE = 0.0
DEFAULT_DIP = 90.0
DEFAULT_RAKE = 0.0
DEFAULT_WIDTH = 0.0
DEFAULT_ZTOR = 0.0


class Origin(object):
    """
    The purpose of this class is to read/store event origin information, which 
    is usually derived from an xml file called "event.xml". There is also a
    mechanism to overwrite information in event.xml with a file in the input
    directory called "source.txt". 

    Note that this class is generally contained within a Rupture instance.

    Event values are: 

        - id: the event id. 
        - created: file creation time (Unix epoch -- seconds since Jan 1, 1970). 
        - lat: hypocenter latitude (decimal degrees; -90 to 90).
        - lon: hypocenter longitude (decimal degrees; -180 to 180).
        - depth: hypocenter depth (km, positive down).
        - locstring: a free-form descriptive string of location.
        - mag: earthquake magnitude.
        - year: 4 digit format.
        - month: 1-12.
        - day: 1-31.
        - hour: 0-23.
        - minute: 0-59.
        - second: 0-59.
        - timezone: abbreviation (i.e., GMT, PST, PDT).
        - mech: a string specifying the rupture mechanism; the accepted types
          are RS, SS, NM, and ALL, for reverse slip, strike slip, normal, and 
          unspecified ruptures, respectively.
    
    For backward-compatibility, we also check for 'type'. If both 'mech' and 
    'type' are missing (or empty strings) then 'mech' is set to ALL.

    """

    def __init__(self, event):
        """
        Construct a Origin object.

        Args:
            event (dict): Dictionary of values. See list above for required keys.

        Returns:
            Origin object.

        """

        #-----------------------------------------------------------------------
        # Check for missing keys
        #-----------------------------------------------------------------------
        missing = []
        for req in REQUIRED_KEYS:
            if req not in list(event.keys()):
                missing.append(req)

        if len(missing):
            raise Exception('Input event dictionary is missing the following '
                            'required keys: "%s"' % (','.join(missing)))

        #-----------------------------------------------------------------------
        # Check some types, ranges, and defaults
        #-----------------------------------------------------------------------
        if not type(event['id']) is str:
            raise Exception('id must be a string.')
        event['lat'] = float(event['lat'])
        if (event['lat'] > 90) or (event['lat'] < -90):
            raise Exception('lat must be between -90 and 90 degrees.')
        event['lon'] = float(event['lon'])
        if (event['lon'] > 180) or (event['lon'] < -180):
            raise Exception('lat must be between -180 and 180 degrees.')
        event['depth'] = float(event['depth'])
        if not type(event['locstring']) is str:
            raise Exception('locstring must be a string.')
        event['mag'] = float(event['mag'])
        if 'mech' in event.keys():
            if event['mech'] == '':
                event['mech'] = DEFAULT_MECH
            if not event['mech'] in RAKEDICT.keys():
                raise Exception('mech must be SS, NM, RS, or ALL.')
        elif 'type' in event.keys():
            event['mech'] = event['type']
            if event['mech'] == '':
                event['mech'] = DEFAULT_MECH
            if not event['mech'] in RAKEDICT.keys():
                raise Exception('mech must be SS, NM, RS, or ALL.')
        else:
            event['mech'] = DEFAULT_MECH

        #-----------------------------------------------------------------------
        # Special handling of time
        #-----------------------------------------------------------------------
        if ('year' in event.keys()) and \
           ('month' in event.keys()) and \
           ('day' in event.keys()) and \
           ('hour' in event.keys()) and \
           ('second' in event.keys()) and \
           ('msec' in event.keys()):
            year = int(event['year'])
            month = int(event['month'])
            day = int(event['day'])
            hour = int(event['hour'])
            minute = int(event['minute'])
            second = float(event['second'])
            msec = int((second - int(second)) * 1e6) # microsec
            second = int(second)
            event['time'] = ShakeDateTime(
                year, month, day, hour, minute, second, msec)
        else:
            event['time'] = ShakeDateTime.utcfromtimestamp(int(time.time()))

        #-----------------------------------------------------------------------
        # Add keys as class attributes
        #-----------------------------------------------------------------------
        for k, v in event.items():
            setattr(self, k, v)

        # What about rake?
        if not hasattr(self, 'rake'):
            if hasattr(self, 'mech'):
                mech = self.mech
                self.rake = RAKEDICT[mech]
            else:
                self.rake = RAKEDICT['ALL']



    @classmethod
    def fromFile(cls, eventxmlfile, sourcefile=None):
        """
        Class method to create a Origin object by specifying an event.xml file,
        a rupture file, and a source.txt file.

        Args:
            eventxmlfile (str): Event xml file (see read_event_file()).
            sourcefile (str): source.txt file (see read_source()).

        Returns:
            Origin object. 

        """
        event = read_event_file(eventxmlfile)
        params = None
        if sourcefile is not None:
            params = read_source(sourcefile)
            event.update(params)
        return cls(event)



    def getHypo(self):
        """
        Returns:
           Hypocenter as OpenQuake Point instance.
        """
        return point.Point(self.lon, self.lat, self.depth)

    def setTectonicRegion(self, region):
        """
        Set tectonic region.

        Args:
            region:
                `TRT <http://docs.openquake.org/oq-hazardlib/master/const.html#openquake.hazardlib.const.TRT>`__ 
                object. 
        """
        if not isinstance(region, TRT):
            raise ValueError(
                'Input tectonic region must be of type openquake.hazardlib.const.TRT')
        self._tectonic_region = copy.deepcopy(region)

    def getTectonicRegion(self):
        """
        Get tectonic region.

        Returns:
            `TRT <http://docs.openquake.org/oq-hazardlib/master/const.html#openquake.hazardlib.const.TRT>`__
            object.
        """
        return copy.deepcopy(self.TectonicRegion)

    def setMechanism(self, mech, rake=None, dip=None):
        """
        Set the earthquake mechanism manually (overriding any values read 
        in from event.xml or source.txt. If rake and dip are not specified, 
        they will be assigned by mechanism as follows:

        +-------+--------+-----+
        | Mech  |  Rake  | Dip |
        +=======+========+=====+
        | RS    |    90  |  40 |
        +-------+--------+-----+
        | NM    |   -90  |  50 |
        +-------+--------+-----+
        | SS    |     0  |  90 |
        +-------+--------+-----+
        | ALL   |    45  |  90 |
        +-------+--------+-----+

        Args:
            mech (str): One of 'RS' (reverse), 'NM' (normal), 'SS' (strike slip), 
                or 'ALL' (unknown).
            rake (float): Value between -360 and 360 degrees. If set, will 
                override default value for mechanism (see table above).
            dip (float): Value betweeen 0 and 90 degrees. If set, will override
                default value for mechanism (see table above). Value will be
                converted to range between -180 and 180 degrees.
        """

        mechs = {'RS': {'rake': 90.0, 'dip': 40.0},
                 'NM': {'rake': -90.0, 'dip': 50.0},
                 'SS': {'rake': 0.0, 'dip': 90.0},
                 'ALL': {'rake': 45.0, 'dip': 90.0}}

        if mech not in list(mechs.keys()):
            raise ShakeMapException(
                'Mechanism must be one of: %s' % str(
                    list(mechs.keys())))

        if dip is not None:
            if dip < 0 or dip > 90:
                raise Exception('Dip must be in range 0-90 degrees.')
            if dip < 0 or dip > 90:
                raise Exception('Dip must be in range 0-90 degrees.')
        else:
            dip = mechs[mech]['dip']

        if rake is not None:
            if rake < -180:
                rake += 360
            if rake > 180:
                rake -= 360
            if rake < -180 or rake > 180:
                raise Exception(
                    'Rake must be transformable to be in range -180 to 180 degrees.')
        else:
            rake = mechs[mech]['rake']

        self.dip = dip
        self.rake = rake
        self.mech = mech


    def computeRepi(self, lon, lat, depth):
        """
        Compute epicentral distance (km). 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            array: Epicentral distance (km).

        """
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        repi = geodetic.distance(self.lon, self.lat, 0.0,
                                 lon, lat, depth)
        repi = repi.reshape(oldshape)
        return repi

    def computeRhyp(self, lon, lat, depth):
        """
        Compute hypocentral distance (km). 

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            array: Hypocentral distance (km).

        """
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        rhyp = geodetic.distance(self.lon, self.lat, self.depth,
                                 lon, lat, depth)
        rhyp = rhyp.reshape(oldshape)
        return rhyp

def read_event_file(eventxml):
    """
    Read event.xml file from disk, returning a dictionary of attributes.
    Input XML format looks like this:

    .. code-block:: xml

         <earthquake 
             id="2008ryan "
             lat="30.9858" 
             lon="103.3639" 
             mag="7.9" 
             year="2008" 
             month="05" 
             day="12" 
             hour="06" 
             minute="28" 
             second="01" 
             timezone="GMT" 
             depth="19.0" 
             locstring="EASTERN SICHUAN, CHINA" 
             created="1211173621" 
             otime="1210573681" 
             type="" 
         />

    Args:
        eventxml (str): Path to event XML file OR file-like object. 
    
    Returns:
       Dictionary with keys as indicated above in earthquake element attributes.

    """

    # fill in default values for mechanism, rake and dip
    # these may be overriden by values in event.xml, source.txt, or by values
    # passed in after reading input data.
#    event = {'mech': DEFAULT_MECH,
#             'rake': DEFAULT_RAKE,
#             'dip': DEFAULT_DIP}

    if isinstance(eventxml, str):
        root = minidom.parse(eventxml)
    else:
        data = eventxml.read()
        root = minidom.parseString(data)

    # Turn XML content into dictionary
    eq = root.getElementsByTagName('earthquake')[0]
    xmldict = dict(eq.attributes.items())
    root.unlink()

    return xmldict



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


