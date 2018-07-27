#!/usr/bin/env python3


# Depth tolerance in km (for determining if top and bottom edges are
# horizontal)
DEPTH_TOL = 0.05

# Maximum ratio of distance off of the plane (relative to edge length) for the
# 4th point to be before being considered non-co-planar and adjusted to
# actually be on the plane?
OFFPLANE_TOLERANCE = 0.05

RAKEDICT = {'SS': 0.0, 'NM': -90.0, 'RS': 90.0, 'ALL': None}

DEFAULT_MECH = 'ALL'
DEFAULT_STRIKE = 0.0
DEFAULT_DIP = 90.0
DEFAULT_RAKE = 0.0
DEFAULT_WIDTH = 0.0
DEFAULT_ZTOR = 0.0

ORIGIN_REQUIRED_KEYS = ['id', 'netid', 'network', 'lat', 'lon', 'depth',
                        'locstring', 'mag', 'time']
# Times can have either integer or floating point (preferred) seconds
TIMEFMT = '%Y-%m-%dT%H:%M:%S.%fZ'
ALT_TIMEFMT = '%Y-%m-%dT%H:%M:%SZ'
