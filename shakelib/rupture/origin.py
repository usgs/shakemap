#!/usr/bin/env python

# stdlib imports
from xml.dom import minidom
import time as time

# third party
from openquake.hazardlib.geo import point

from impactutils.time.ancient_time import HistoricTime
from shakelib.utils.exception import ShakeLibException
from shakelib.rupture import constants

TIMEFMT = '%Y-%m-%dT%H:%M:%SZ'

class Origin(object):
    """
    The purpose of this class is to read/store event origin information, which
    is usually derived from an xml file called "event.xml". There is also a
    mechanism to overwrite information in event.xml with a file in the input
    directory called "source.txt".

    Note that this class is generally contained within a Rupture instance.

    Event values are:

        - eventsourcecode: the event id.
        - created: file creation time (Unix epoch - seconds since Jan 1, 1970).
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
        Construct an Origin object.

        Args:
            event (dict): Dictionary of values. See list above for required
            keys.

        Returns:
            Origin object.
        Raises:
            ValueError: When input time is not and cannot be converted to
                HistoricTime object.
        """

        # ---------------------------------------------------------------------
        # Check for missing keys
        # ---------------------------------------------------------------------
        missing = []
        for req in constants.ORIGIN_REQUIRED_KEYS:
            if req not in list(event.keys()):
                missing.append(req)

        if len(missing):
            raise Exception('Input event dictionary is missing the following '
                            'required keys: "%s"' % (','.join(missing)))

        # ---------------------------------------------------------------------
        # Check some types, ranges, and defaults
        # ---------------------------------------------------------------------
        if not type(event['eventsourcecode']) is str:
            raise Exception('eventsourcecode must be a string.')

        if (event['lat'] > 90) or (event['lat'] < -90):
            raise Exception('lat must be between -90 and 90 degrees.')

        if (event['lon'] > 180) or (event['lon'] < -180):
            raise Exception('lat must be between -180 and 180 degrees.')

        # make sure that time is an HistoricTime instance
        if 'time' in event:
            if isinstance(event['time'],str):
                try:
                    event['time'] = HistoricTime.strptime(event['time'],TIMEFMT)
                except ValueError:
                    fmt = 'Input time string %s cannot be converted to datetime.'
                    raise ValueError(fmt % event['time'])

        if 'mech' in event.keys():
            if event['mech'] == '':
                event['mech'] = constants.DEFAULT_MECH
            if not event['mech'] in constants.RAKEDICT.keys():
                raise Exception('mech must be SS, NM, RS, or ALL.')
        elif 'type' in event.keys():
            event['mech'] = event['type']
            if event['mech'] == '':
                event['mech'] = constants.DEFAULT_MECH
            if not event['mech'] in constants.RAKEDICT.keys():
                raise Exception('mech must be SS, NM, RS, or ALL.')
        else:
            event['mech'] = constants.DEFAULT_MECH


        # ---------------------------------------------------------------------
        # Add keys as class attributes
        # ---------------------------------------------------------------------
        for k, v in event.items():
            if k == 'rake':
                setattr(self, k, float(v))
            else:
                setattr(self, k, v)

        # What about rake?
        if not hasattr(self, 'rake'):
            if hasattr(self, 'mech'):
                mech = self.mech
                self.rake = constants.RAKEDICT[mech]
            else:
                self.rake = constants.RAKEDICT['ALL']

        if self.rake is None:
            self.rake = 0.0

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

    def setMechanism(self, mech=None, rake=None, dip=None):
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
            mech (str): One of 'RS' (reverse), 'NM' (normal), 'SS' (strike
                slip), or 'ALL' (unknown).
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
            raise ShakeLibException(
                'Mechanism must be one of: %s' % str(
                    list(mechs.keys())))

        if dip is not None:
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
                    'Rake must be transformable to be in range -180 to 180 '
                    'degrees.')
        else:
            rake = mechs[mech]['rake']

        self.dip = dip
        self.rake = rake
        self.mech = mech


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
       dict: Dictionary with keys:
         - eventsourcecode Origin network and origin code (i.e., us2017abcd).
         - eventsource Origin network ("us").
         - time Origin time as an HistoricTime object.
         - lat Origin latitude
         - lon Origin longitude
         - depth Origin depth
         - mag Origin magnitude
         - created Process time as an HistoricTime object.
         - locstring Location string
         - mechanism Moment mechanism, one of:
           - 'RS' (Reverse)
           - 'SS' (Strike-Slip)
           - 'NM' (Normal)
           - 'ALL' (Undetermined)
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

    eqdict = {}
    eqdict['eventsourcecode'] = xmldict['id']
    if 'network' in xmldict:
        eqdict['eventsource'] = xmldict['network']
    else:
        eqdict['eventsource'] = 'us' #??

    #look for the productcode attribute
    if 'productcode' in xmldict:
        eqdict['productcode'] = xmldict['productcode']

    # fix eventsourcecode if not specified correctly
    if not eqdict['eventsourcecode'].startswith(eqdict['eventsource']):
        eqdict['eventsourcecode'] = eqdict['eventsource'] + eqdict['eventsourcecode']

    year = int(xmldict['year'])
    month = int(xmldict['month'])
    day = int(xmldict['day'])
    hour = int(xmldict['hour'])
    minute = int(xmldict['minute'])
    second = int(xmldict['second'])
    microseconds = int((second - int(xmldict['second']))*1e6)
    eqdict['time'] = HistoricTime(year,month,day,hour,minute,second,microseconds)
    eqdict['lat'] = float(xmldict['lat'])
    eqdict['lon'] = float(xmldict['lon'])
    eqdict['depth'] = float(xmldict['depth'])
    eqdict['mag'] = float(xmldict['mag'])

    # make created field in event.xml optional - set to current UTC time if not
    # supplied.
    if 'created' in xmldict:
        eqdict['created'] = HistoricTime.utcfromtimestamp(int(xmldict['created']))
    else:
        eqdict['created'] = HistoricTime.utcnow()

    eqdict['locstring'] = xmldict['locstring']
    
    if 'mech' in xmldict:
        eqdict['mech'] = xmldict['mech']
    return eqdict


def read_source(sourcefile):
    """
    Read source.txt file, which has lines like key=value.

    Args:
        sourcefile: Path to source.txt file OR file-like object

    Returns:
        dict: Dictionary containing key/value pairs from file.
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
