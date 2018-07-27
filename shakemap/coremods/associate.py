# stdlib imports
import os.path

# third party imports

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths
from shakemap.utils.amps import AmplitudeHandler
from shakelib.rupture.origin import read_event_file
from shakelib.rupture import constants


class AssociateModule(CoreModule):
    """
    associate -- Associate amps in the database with the event, and write
                 XML data file to the event's current directory.
    """

    command_name = 'associate'

    def execute(self):
        """
        Associate amps and write unassoc_<datetime>_dat.xml.
        """
        install_path, data_path = get_config_paths()

        amp_handler = AmplitudeHandler(install_path, data_path)

        event = amp_handler.getEvent(self._eventid)
        if event is None:
            #
            # This shouldn't ever happen, but the code is here just
            # in case it does
            #
            datadir = os.path.join(data_path, self._eventid, 'current')
            if not os.path.isdir(datadir):
                raise NotADirectoryError('%s is not a valid directory.' %
                                         datadir)
            eventxml = os.path.join(datadir, 'event.xml')
            if not os.path.isfile(eventxml):
                raise FileNotFoundError('%s does not exist.' % eventxml)
            origin = read_event_file(eventxml)

            event = {'id': self._eventid,
                     'netid': origin['netid'],
                     'network': origin['network'],
                     'time': origin['time'].strftime(constants.TIMEFMT),
                     'lat': origin['lat'],
                     'lon': origin['lon'],
                     'depth': origin['depth'],
                     'mag': origin['mag'],
                     'locstring': origin['locstring']}
            amp_handler.insertEvent(event)

        amp_handler.associateOne(self._eventid, pretty_print=True)
