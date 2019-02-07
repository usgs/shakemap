# stdlib imports
import os.path
import glob
import logging

# third party imports

# local imports
from .base import CoreModule
from shakemap.utils.logging import get_logging_config
from impactutils.io.smcontainers import ShakeMapOutputContainer
from shakemap.utils.config import get_config_paths


class HistoryModule(CoreModule):
    """
    history -- Output the version history of an event.
    """

    command_name = 'history'

    def __init__(self, eventid):
        """
        Instantiate a ContourModule class with an event ID.
        """
        self._eventid = eventid
        log_config = get_logging_config()
        log_name = log_config['loggers'].keys()[0]
        self.logger = logging.getLogger(log_name)

    def execute(self):
        """
        Output the version history of an event.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
        """
        _, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        backups = glob.glob(os.path.join(data_path, self._eventid, 'backup*'))
        backups.sort(reverse=True)
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)

        # First try the current results file...
        datafile = os.path.join(datadir, 'products', 'shake_result.hdf')
        if os.path.isfile(datafile):
            # Open the ShakeMapOutputContainer and extract the data
            container = ShakeMapOutputContainer.load(datafile)
            try:
                metadata = container.getMetadata()
            except LookupError:
                print("\nNo version history available for this event.\n")
                return
            history = (metadata['processing']['shakemap_versions']
                               ['map_data_history'])
            final = False
            if len(backups) > 0:
                last_ver = int(backups[0][-4:])
                last_hist = history[-1][2]
                if last_ver == last_hist:
                    final = True
            print_history(history, final=final)
            return

        # Nope. Are there any backup files?
        if len(backups) == 0:
            print("\nNo version history available for this event.\n")
            return

        # There should be a results file in the backup directory...
        datafile = os.path.join(data_path, self._eventid, backups[0],
                                'products', 'shake_result.hdf')
        if os.path.isfile(datafile):
            # Open the ShakeMapOutputContainer and extract the data
            container = ShakeMapOutputContainer.load(datafile)
            try:
                metadata = container.getMetadata()
            except LookupError:
                print("\nNo version history available for this event.\n")
                return
            history = (metadata['processing']['shakemap_versions']
                               ['map_data_history'])
            print_history(history, final=True)
            return

        print("\nNo version history available for this event.\n")
        return


def print_history(history, final=False):

    if len(history) == 0:
        print("\nVersion history is empty.\n")
        return

    print("\nVersion history:")
    print("Timestamp | Originator | Version | Comment")
    print("------------------------------------------")
    for ix, line in enumerate(history):
        if final is False and ix == len(history) - 1:
            asterisk = '*'
        else:
            asterisk = ''
        print("%s | %s | %d | %s%s" % (line[0], line[1], line[2], line[3],
                                       asterisk))
    print("")
    if final is False:
        print("*Not finalized.\n")
