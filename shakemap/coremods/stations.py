# stdlib imports
import os.path
import json

# third party imports
from shakelib.utils.containers import ShakeMapOutputContainer
from configobj import ConfigObj

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths

ALLLOWED_FORMATS = ['json']


class StationModule(CoreModule):
    """
    stations -- Generate stationlist.json from shake_result.hdf.
    """

    command_name = 'stations'

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

        # get the path to the products.conf file, load the config
        config_file = os.path.join(install_path, 'config', 'products.conf')
        config = ConfigObj(config_file)

        # create ShakeMap station data file
        formats = config['products']['info']['formats']
        for fformat in formats:
            if fformat not in ALLLOWED_FORMATS:
                self.logger.warn('Specified format %s not in list of defined '
                                 'formats.  Skipping.' % fformat)
                continue
            if fformat == 'json':
                self.logger.debug('Writing rupture.json file...')
                station_dict = container.getStationDict()
                station_file = os.path.join(datadir,'stationlist.json')
                f = open(station_file,'w')
                json.dump(station_dict, f)
                f.close()
