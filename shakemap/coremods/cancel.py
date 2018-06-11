# stdlib imports
import os.path
from datetime import datetime

# third party imports
from shakelib.utils.containers import ShakeMapOutputContainer
from configobj import ConfigObj


# local imports
from .base import CoreModule
from .transfer import _transfer
from shakemap.utils.config import get_config_paths, get_data_path
from shakemap.utils.amps import AmplitudeHandler

TIMEFMT = '%Y-%m-%dT%H:%M:%SZ'


class CancelModule(CoreModule):
    """
    cancel -- Send a cancel message via any configured transfer mechanisms.
    """

    command_name = 'cancel'

    def execute(self):
        """
        Cancel ShakeMap products using methods configured in transfer.conf.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        products_dir = os.path.join(datadir, 'products')
        if not os.path.isdir(products_dir):
            raise NotADirectoryError('%s does not exist.' % products_dir)

        # get the path to the transfer.conf spec file
        configspec = os.path.join(get_data_path(), 'transferspec.conf')

        # look for an event specific transfer.conf file
        transfer_conf = os.path.join(datadir, 'transfer.conf')
        if not os.path.isfile(transfer_conf):
            # if not there, use the system one
            transfer_conf = os.path.join(
                install_path, 'config', 'transfer.conf')
            if not os.path.isfile(transfer_conf):
                raise FileNotFoundError('%s does not exist.' % transfer_conf)

        # get the config information for transfer
        config = ConfigObj(transfer_conf, configspec=configspec)

        # get the output container with all the things in it
        datafile = os.path.join(products_dir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the ShakeMapOutputContainer and extract the data
        container = ShakeMapOutputContainer.load(datafile)

        # call the transfer method
        self.logger.info('Sending cancel message...')
        _transfer(config, container, products_dir, cancel=True)

        # Create a file called CANCEL in the data directory. The
        # shake program will look for this and not run if present.
        self.logger.info('Creating cancel file...')
        cancelfile = os.path.join(datadir, 'CANCEL')
        with open(cancelfile, 'wt') as cfile:
            cfile.write('Event cancelled at %s\n' %
                        datetime.utcnow().strftime(TIMEFMT))

        # delete the event from the database
        handler = AmplitudeHandler(install_path, data_path)
        handler.deleteEvent(self._eventid)
        container.close()
