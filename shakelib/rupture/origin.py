#!/usr/bin/env python

# stdlib imports
import warnings
import os.path

# third party
from lxml import etree
import defusedxml.ElementTree as dET
from openquake.hazardlib.geo import point
from obspy.io.quakeml.core import _is_quakeml
from obspy.core.event import read_events
from strec.tensor import fill_tensor_from_components

from impactutils.time.ancient_time import HistoricTime
from shakelib.utils.exception import ShakeLibException
from shakelib.rupture import constants


class Origin(object):
    """The purpose of this class is to read/store event origin information,
    which is usually derived from an xml file called "event.xml". There is
    also a mechanism to overwrite information in event.xml with a file in
    the input directory called "source.txt".

    Note that this class is generally contained within a Rupture instance.

    Event values are:

        - id: the event id.
        - netid: the source network code
        - network: a string describing the source network
        - lat: hypocenter latitude (decimal degrees; -90 to 90).
        - lon: hypocenter longitude (decimal degrees; -180 to 180).
        - depth: hypocenter depth (km, positive down).
        - locstring: a free-form descriptive string of location.
        - mag: earthquake magnitude.
        - time: a HistoricTime object containing the origin date/time
        - mech: a string specifying the rupture mechanism; the accepted types
          are RS, SS, NM, and ALL, for reverse slip, strike slip, normal, and
          unspecified ruptures, respectively.
        - reference: A string describing the data source (optional)
        - productcode: The PDL product code for this ShakeMap

    For backward-compatibility, we also check for 'type'. If both 'mech' and
    'type' are missing (or empty strings) then 'mech' is set to ALL.
    """

    def __init__(self, event):
        """ Construct an Origin object.

        Args:
            event (dict): Dictionary of values. See list above for required
            keys.

        Returns:
            Origin object.

        Raises:
            KeyError: When one or more required values are missing from the
                      event object.
            TypeError: If the event ID is not a string.
            ValueError: If the latitude or longitude are outside of the
                        allowed ranges, or if the 'mech' value is not
                        an allowed type.
        """

        # ---------------------------------------------------------------------
        # Check for missing keys
        # ---------------------------------------------------------------------
        missing = []
        ekeys = list(event.keys())
        for req in constants.ORIGIN_REQUIRED_KEYS:
            if req not in ekeys:
                missing.append(req)

        if len(missing):
            raise KeyError('Input event dictionary is missing the following '
                           'required keys: "%s"' % (', '.join(missing)))

        # ---------------------------------------------------------------------
        # Check some types, ranges, and defaults
        # ---------------------------------------------------------------------
        if not type(event['id']) is str and event['id'] is not None:
            raise TypeError('event id  must be a string.')

        if float(event['lat']) > 90 or float(event['lat']) < -90:
            raise ValueError('lat must be between -90 and 90 degrees.')

        if float(event['lon']) > 180 or float(event['lon']) < -180:
            raise ValueError('lon must be between -180 and 180 degrees.')

        if 'type' in ekeys and 'mech' not in ekeys:
            event['mech'] = event['type']
            del event['type']
        if 'mech' in event.keys():
            if event['mech'] == '':
                event['mech'] = constants.DEFAULT_MECH
            if not event['mech'] in constants.RAKEDICT.keys():
                raise ValueError('mech must be SS, NM, RS, or ALL.')
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
            self.rake = constants.RAKEDICT[self.mech]

        if self.rake is None:
            self.rake = 0.0

        return

    @classmethod
    def fromFile(cls, eventxmlfile, sourcefile=None, momentfile=None):
        """Class method to create a Origin object by specifying an event.xml
        file, a rupture file, a source.txt file, and a moment.xml file.

        Args:
            eventxmlfile (str): Event xml file (see read_event_file()).
            sourcefile (str): source.txt file (see read_source()).
            momentfile (str): moment.xml file (see read_moment_quakeml()).

        Returns:
            Origin object.
        """
        event = read_event_file(eventxmlfile)
        params = None
        if sourcefile is not None:
            params = read_source(sourcefile)
            event.update(params)
        if momentfile is not None:
            params = read_moment_quakeml(momentfile)
            event.update(params)
        return cls(event)

    def getHypo(self):
        """
        Returns:
           Hypocenter as OpenQuake Point instance.
        """
        return point.Point(self.lon, self.lat, self.depth)

    def setMechanism(self, mech=None, rake=None, dip=None):
        """Set the earthquake mechanism manually (overriding any values read
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
        return


def write_event_file(event, xmlfile):
    """Write event.xml file.

    Args:
       event (dict): Dictionary containing fields [field, (type, example)]:

           - id (str, "us2008abcd")
           - netid (str, "us")
           - network (str, "USGS National Network")
           - lat (float, 42.1234)
           - lon (float, -85.1234)
           - depth (float, 24.1)
           - mag (float, 7.9)
           - time (HistoricTime object)
           - locstring (str, "East of the Poconos")
           - mech (str, "RS", "SS", "NM", or 'ALL') (optional)
           - reference (str, data source: 'Smith et al. 2016') (optional)
           - event_type (str, 'ACTUAL' or 'SCENARIO') (optional)
           - productcode (str, 'us2000wxyz_zoom')

    Returns:
        nothing: Nothing.
    """
    try:
        root = etree.Element('earthquake')
        root.attrib['id'] = event['id']
        root.attrib['netid'] = event['netid']
        root.attrib['network'] = event['network']
        root.attrib['lat'] = '%.4f' % event['lat']
        root.attrib['lon'] = '%.4f' % event['lon']
        root.attrib['depth'] = '%.1f' % event['depth']
        root.attrib['mag'] = '%.1f' % event['mag']
        root.attrib['time'] = event['time'].strftime(constants.ALT_TIMEFMT)
        root.attrib['locstring'] = event['locstring']
        if 'mech' in event:
            root.attrib['mech'] = event['mech']
        if 'reference' in event:
            root.attrib['reference'] = event['reference']
        if 'productcode' in event:
            root.attrib['productcode'] = event['productcode']
        if 'event_type' in event:
            if event['event_type'] not in ['ACTUAL', 'SCENARIO']:
                raise AttributeError(
                    'event_type must be "ACTUAL" or "SCENARIO"')
            root.attrib['event_type'] = event['event_type']
        else:
            root.attrib['event_type'] = 'ACTUAL'

        tree = etree.ElementTree(root)
        tree.write(xmlfile, pretty_print=True)
    except Exception as e:
        return str(e)

    return None


def read_event_file(eventxml):
    """Read event.xml file from disk, returning a dictionary of attributes.
    Input XML format looks like this (all elements are required unless
    explicitly labeled optional):

    .. code-block:: xml

         <earthquake
             id="2008ryan"
             netid="us"
             network="USGS National Network" (required but may be empty string)
             lat="30.9858"
             lon="103.3639"
             mag="7.9"
             depth="19.0"
             time="YYYY-mm-ddTHH:MM:SS.ffffffZ" (omitting fractional seconds
                                                is also supported)
             locstring="EASTERN SICHUAN, CHINA"
             mech='SS' | 'NM' | 'RS' | 'ALL' (optional)
             reference="Smith et al. 2016" (optional)
             productcode='us2008ryan' (optional)
         />

    Args:
        eventxml (str): Path to event XML file OR file-like object.

    Returns:
       dict: Dictionary with keys:
         - id: Origin network and origin code (i.e., us2017abcd).
         - netid: Origin network ("us").
         - network: (A long-form description of the network)
         - lat: Origin latitude
         - lon: Origin longitude
         - mag: Origin magnitude
         - depth: Origin depth
         - time: Origin time as an HistoricTime object.
         - locstring: Location string
         - mech: (optional) Moment mechanism, one of:
           - 'RS' (Reverse)
           - 'SS' (Strike-Slip)
           - 'NM' (Normal)
           - 'ALL' (Undetermined)
         - reference: (optional) A description of the source of the data.
         - productcode: (optional) This product source's unique code for this
                        particular ShakeMap.

    Raises:
        ValueError: If the time string cannot be parsed into a datetime object
        KeyError: If any of the required attributes are missing from event.xml
    """
    if isinstance(eventxml, str):
        tree = dET.parse(eventxml)
        root = tree.getroot()
    else:
        data = eventxml.read()
        root = dET.fromstring(data)

    # Turn XML content into dictionary
    if root.tag == 'earthquake':
        xmldict = dict(root.items())
    else:
        eq = root.find('earthquake')
        xmldict = dict(eq.items())

    eqdict = {}

    #########################################################
    # A Short Primer on PDL-style Identifiers
    # Because Everybody (Including The Author) Forgets It.
    #
    # In PDL, there are 4 identifiers that fully specify a product:
    # - source The network that generated the *product* (us, ci, etc.).
    # - code   The unique ID string that identifies this product,
    #          usually prepended by *source* (us2008abcd).
    # - eventsource The network that created the *origin* (us, ci, etc.)
    # - eventsourcecode The code within that network that uniquely
    #                   identifies the event (2008abcd).
    #
    # For our purposes, we're storing *source* and *code* as
    # *productsource* and *productcode* respectively in the
    # container, in an effort to reduce confusion about their
    # meaning. Time will tell.
    #########################################################

    # read in the id fields
    eqdict['id'] = xmldict['id']
    # This isn't optional, but maybe it isn't in some old files
    if 'network' in xmldict:
        eqdict['network'] = xmldict['network']
    else:
        eqdict['network'] = ""

    eqdict['netid'] = xmldict['netid']

    # look for the productcode attribute in the xml,
    # otherwise use the event id
    if 'productcode' in xmldict:
        eqdict['productcode'] = xmldict['productcode']
    elif isinstance(eventxml, str):
        eqdict['productcode'] = eqdict['id']
    else:
        # It's up to the user of this data how to construct the
        # product code
        pass

    # Support old event file date/times
    if 'time' in xmldict:
        try:
            eqdict['time'] = HistoricTime.strptime(
                xmldict['time'],
                constants.TIMEFMT)
        except ValueError:
            try:
                eqdict['time'] = HistoricTime.strptime(
                    xmldict['time'],
                    constants.ALT_TIMEFMT)
            except ValueError:
                raise ValueError("Couldn't convert %s to HistoricTime" %
                                 xmldict['time'])
    else:
        if 'year' not in xmldict or 'month' not in xmldict or \
           'day' not in xmldict or 'hour' not in xmldict or \
           'minute' not in xmldict or 'second' not in xmldict:
            raise ValueError("Missing date/time elements in event file.")
        eqdict['time'] = HistoricTime.datetime(
            xmldict['year'],
            xmldict['month'],
            xmldict['day'],
            xmldict['hour'],
            xmldict['minute'],
            xmldict['second'])

    eqdict['lat'] = float(xmldict['lat'])
    eqdict['lon'] = float(xmldict['lon'])
    eqdict['depth'] = float(xmldict['depth'])
    eqdict['mag'] = float(xmldict['mag'])
    eqdict['locstring'] = xmldict['locstring']

    if 'mech' in xmldict:
        eqdict['mech'] = xmldict['mech']
    # Older files may have "type" instead of "mech"
    if 'type' in xmldict:
        eqdict['type'] = xmldict['type']
    if 'reference' in xmldict:
        eqdict['reference'] = xmldict['reference']

    return eqdict


def read_moment_quakeml(momentfile):
    """Read moment parameters from a QuakeML file.

    Args:
        momentfile (str): Path to a QuakeML file.

    Returns:
        dict: Empty if momentfile is somehow not valid, or:
              - moment:
                  - T:
                    - azimuth
                    - plunge
                  - N:
                    - azimuth
                    - plunge
                  - P:
                    - azimuth
                    - plunge
                  - NP1:
                    - strike: float
                    - dip: float
                    - rake: float
                  - NP2:
                    - strike: float
                    - dip: float
                    - rake: float
                  - type: string
                  - source: string
    """
    params = {}
    if not _is_quakeml(momentfile):
        return params

    # obspy spits out a bunch of warnings we don't care about
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        catalog = read_events(momentfile)

    if not len(catalog.events):
        return params
    event = catalog.events[0]
    if not len(event.focal_mechanisms):
        return params
    focal = catalog.events[0].focal_mechanisms[0]

    # get the source and type information, if possible
    mtype = 'unknown'
    if focal.method_id is not None:
        mtype = focal.method_id

    msource = 'unknown'
    if focal.creation_info is not None:
        if focal.creation_info.agency_id is not None:
            msource = focal.creation_info.agency_id

    if focal.nodal_planes is None:
        return params
    if focal.nodal_planes.nodal_plane_1 is None:
        if focal.moment_tensor.tensor.m_rr is not None:
            mrr = focal.moment_tensor.tensor.m_rr
            mtt = focal.moment_tensor.tensor.m_tt
            mpp = focal.moment_tensor.tensor.m_pp
            mrt = focal.moment_tensor.tensor.m_rt
            mrp = focal.moment_tensor.tensor.m_rp
            mtp = focal.moment_tensor.tensor.m_tp
            params = fill_tensor_from_components(
                mrr, mtt,
                mpp, mrt,
                mrp, mtp,
                source=msource,
                mtype=mtype)
    else:
        plane1 = focal.nodal_planes.nodal_plane_1
        plane2 = focal.nodal_planes.nodal_plane_2
        params['NP1'] = {
            'strike': plane1.strike,
            'dip': plane1.dip,
            'rake': plane1.rake
        }
        params['NP2'] = {
            'strike': plane2.strike,
            'dip': plane2.dip,
            'rake': plane2.rake
        }
        params['T'] = {
            'azimuth': focal.principal_axes.t_axis.azimuth,
            'plunge': focal.principal_axes.t_axis.plunge
        }
        params['N'] = {
            'azimuth': focal.principal_axes.n_axis.azimuth,
            'plunge': focal.principal_axes.n_axis.plunge
        }
        params['P'] = {
            'azimuth': focal.principal_axes.p_axis.azimuth,
            'plunge': focal.principal_axes.p_axis.plunge
        }

    moment = {'moment': params}
    return moment


def read_source(sourcefile):
    """Read source.txt file, which has lines like key=value.

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
