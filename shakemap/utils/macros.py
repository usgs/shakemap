# from local imports
# from shakemap.utils.queue import TIMEFMT

# stdlib imports
from datetime import datetime

from shakelib.rupture import constants

DATE_STR_FMT = '%b %d, %Y'
TIME_OF_DAY_FMT = '%H:%M:%S'


def get_macros(info):
    """Return a dictionary containing macros that can be used in shakemail
    text.

    Args:
        info (dict): Dictionary returned from
            ShakeMapOutputContainer.getMetadata().
    Returns:
        dict: Dictionary containing following fields:
              - MAG Event magnitude.
              - LOC Location string.
              - LAT Event latitude.
              - LON Event longitude.
              - DEP Event depth.
              - DATETIME Event date/time (i.e., 2018-01-18T11:34:25.123456)
              - DATE Event date (i.e., "Jan 31, 2018")
              - TIME Event time (i.e., "11:34:23")
              - VERSION ShakeMap map version (i.e., 1, 2, 3, etc.)
              - EVENTID Earthquake event ID.
              - PRODUCT_CODE Unique code describing the ShakeMap product.
              - NETID Earthquake network ID.
    """
    macros = {}
    macros['MAG'] = info['input']['event_information']['magnitude']
    macros['LOC'] = info['input']['event_information']['location']
    macros['LAT'] = info['input']['event_information']['latitude']
    macros['LON'] = info['input']['event_information']['longitude']
    macros['DEP'] = info['input']['event_information']['depth']
    macros['DATETIME'] = info['input']['event_information']['origin_time']
    try:
        dtime = datetime.strptime(macros['DATETIME'], constants.TIMEFMT)
    except ValueError:
        dtime = datetime.strptime(macros['DATETIME'], constants.ALT_TIMEFMT)
    macros['DATE'] = dtime.strftime(DATE_STR_FMT)
    macros['TIME'] = dtime.strftime(TIME_OF_DAY_FMT)
    macros['VERSION'] = str(
        info['processing']['shakemap_versions']['map_version'])
    macros['EVENTID'] = info['input']['event_information']['event_id']
    macros['PRODUCT_CODE'] = info['input']['event_information']['productcode']
    if 'netid' in info['input']['event_information']:
        macros['NETID'] = info['input']['event_information']['netid']
    else:
        macros['NETID'] = ''

    return macros
