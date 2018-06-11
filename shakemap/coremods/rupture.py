# stdlib imports
import os.path
import json
from collections import OrderedDict

# third party imports
from shakelib.utils.containers import ShakeMapOutputContainer

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths

ALLOWED_FORMATS = ['json']


class RuptureModule(CoreModule):
    """
    rupture -- Generate rupture.json from shake_result.hdf.
    """

    command_name = 'rupture'
    targets = [r'products/rupture\.json']
    dependencies = [('products/shake_result.hdf', True)]

    # supply here a data structure with information about files that
    # can be created by this module.
    contents = OrderedDict.fromkeys(['ruptureJSON'])
    contents['ruptureJSON'] = {
        'title': 'Fault Rupture',
        'caption': 'JSON Representation of Fault Rupture.',
        'formats': [{'filename': 'rupture.json',
                     'type': 'application/json'}]
    }

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

        # create ShakeMap rupture file
        for fformat in ALLOWED_FORMATS:
            if fformat == 'json':
                self.logger.info('Writing rupture.json file...')
                rupture_dict = container.getRuptureDict()
                rupture_file = os.path.join(datadir, 'rupture.json')
                f = open(rupture_file, 'w')
                json.dump(rupture_dict, f)
                f.close()

        container.close()
