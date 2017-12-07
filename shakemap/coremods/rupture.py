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


class RuptureModule(CoreModule):
    """
    rupture -- Generate rupture.json from shake_result.hdf.
    """

    command_name = 'rupture'

    def execute(self):
        """
        Write rupture.json file.

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

        # create ShakeMap rupture file
        formats = config['products']['info']['formats']
        for fformat in formats:
            if fformat not in ALLLOWED_FORMATS:
                self.logger.warn('Specified format %s not in list of defined '
                                 'formats.  Skipping.' % fformat)
                continue
            if fformat == 'json':
                self.logger.info('Writing rupture.json file...')
                rupture_dict = container.getRuptureDict()
                rupture_file = os.path.join(datadir, 'rupture.json')
                f = open(rupture_file, 'w')
                json.dump(rupture_dict, f)
                f.close()
