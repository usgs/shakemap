# stdlib imports
import os.path
from collections import OrderedDict

# third party imports
from shakelib.utils.containers import ShakeMapOutputContainer

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths


class InfoModule(CoreModule):
    """
    info -- Extract info.json from shake_result.hdf and write it as a file.
    """

    command_name = 'info'

    # supply here a data structure with information about files that
    # can be created by this module.
    contents = OrderedDict.fromkeys(['supplementalInformation'])
    cap = 'ShakeMap processing parameters and map summary information.'
    contents['supplementalInformation'] = {'title':'Supplemental Information',
                             'caption':cap,
                             'formats':[{'filename':'info.json',
                                         'type':'application/json'}
                                       ]
                            }

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

        # create ShakeMap metadata file
        self.logger.debug('Writing info.json file...')
        infostring = container.getString('info.json')
        info_file = os.path.join(datadir, 'info.json')
        f = open(info_file, 'wt')
        f.write(infostring)
        f.close()

