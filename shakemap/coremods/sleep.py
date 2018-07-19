"""
Sleep for a specified number of seconds.
"""

# stdlib imports
import argparse
import inspect
import logging
import time

# third party imports

# local imports
from .base import CoreModule
from shakemap.utils.config import get_logging_config


class SleepModule(CoreModule):
    """
    sleep -- Sleep for a number of seconds.
    """

    command_name = 'sleep'

    def __init__(self, eventid, seconds=None):
        """
        Instantiate a CoreModule class with an event ID.
        """
        self._eventid = eventid
        log_config = get_logging_config()
        log_name = log_config['loggers'].keys()[0]
        self.logger = logging.getLogger(log_name)
        if seconds is not None:
            self.seconds = seconds

    def execute(self):
        """
        Sleep for the specified number of seconds and return. The default
        is 60 seconds.
        """

        # Prompt for a comment string if none is provided on the command line
        if self.seconds is None:
            self.seconds = 60

        time.sleep(self.seconds)

    def parseArgs(self, arglist):
        """
        Set up the object to accept the --seconds flag.
        """
        parser = argparse.ArgumentParser(
            prog=self.__class__.command_name,
            description=inspect.getdoc(self.__class__))
        parser.add_argument('-s', '--seconds', help='Specify the number of '
                            'seconds that the module should sleep.')
        #
        # This line should be in any modules that overrides this
        # one. It will collect up everything after the current
        # modules options in args.rem, which should be returned
        # by this function. Note: doing parser.parse_known_args()
        # will not work as it will suck up any later modules'
        # options that are the same as this one's.
        #
        parser.add_argument('rem', nargs=argparse.REMAINDER,
                            help=argparse.SUPPRESS)
        args = parser.parse_args(arglist)
        self.seconds = int(args.seconds)
        return args.rem
