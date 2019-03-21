# stdlib imports
import os.path
import json

# third party imports
import numpy as np
from impactutils.io.smcontainers import ShakeMapOutputContainer

# local imports
from .base import CoreModule, Contents
from shakemap.utils.config import get_config_paths


class InfoModule(CoreModule):
    """
    info -- Extract info.json from shake_result.hdf and write it as a file.
    """

    command_name = 'info'
    targets = [r'products/info\.json']
    dependencies = [('products/shake_result.hdf', True)]

    def __init__(self, eventid):
        super(InfoModule, self).__init__(eventid)
        self.contents = Contents(None, None, eventid)

    def execute(self):
        """
        Write info.json metadata file.

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

        # Create ShakeMap metadata file
        self.logger.debug('Writing info.json file...')
        info = container.getMetadata()

        # Clean up strec results to be valid json
        if 'strec' in info:
            for k, v in info['strec'].items():
                if isinstance(v, float):
                    if not np.isfinite(v):
                        info['strec'][k] = None

        infostring = json.dumps(info, allow_nan=False)
        info_file = os.path.join(datadir, 'info.json')
        f = open(info_file, 'wt')
        f.write(infostring)
        f.close()
        container.close()
        cap = 'ShakeMap processing parameters and map summary information.'
        self.contents.addFile('supplementalInformation',
                              'Supplemental Information', cap,
                              'info.json', 'application/json')
