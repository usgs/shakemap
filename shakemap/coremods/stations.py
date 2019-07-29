# stdlib imports
import os.path
import json

# third party imports
from impactutils.io.smcontainers import ShakeMapOutputContainer

# local imports
from .base import CoreModule, Contents
from shakemap.utils.config import get_config_paths

ALLOWED_FORMATS = ['json']


class StationModule(CoreModule):
    """
    stations -- Generate stationlist.json from shake_result.hdf.
    """

    command_name = 'stations'
    targets = [r'products/stationlist\.json']
    dependencies = [('products/shake_result.hdf', True)]

    def __init__(self, eventid):
        super(StationModule, self).__init__(eventid)
        self.contents = Contents(None, None, eventid)

    def execute(self):
        """Write stationlist.json file.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        datafile = os.path.join(datadir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the ShakeMapOutputContainer and extract the data
        container = ShakeMapOutputContainer.load(datafile)

        # create ShakeMap station data file
        for fformat in ALLOWED_FORMATS:
            if fformat == 'json':
                self.logger.debug('Writing rupture.json file...')
                station_dict = container.getStationDict()
                station_file = os.path.join(datadir, 'stationlist.json')
                f = open(station_file, 'w')
                json.dump(station_dict, f)
                f.close()

        container.close()

        self.contents.addFile('stationJSON', 'Station List',
                              'List of ShakeMap input data.',
                              'stationlist.json', 'application/json')
