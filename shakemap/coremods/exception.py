"""
Raise an exception (used for testing.)
"""

# stdlib imports
import logging

# third party imports

# local imports
from .base import CoreModule
from shakemap.utils.logging import get_logging_config


class ExceptionModule(CoreModule):
    """
    Module to raise an exception for testing purposes.
    """

    command_name = "exception"

    def __init__(self, eventid, seconds=None):
        """
        Instantiate a CoreModule class with an event ID.
        """
        self._eventid = eventid
        log_config = get_logging_config()
        log_name = log_config["loggers"].keys()[0]
        self.logger = logging.getLogger(log_name)
        if seconds is not None:
            self.seconds = seconds

    def execute(self):
        """
        Raise an Exception object.

        This module exists for the purposes of testing shake's exception
        handling logic.
        """

        raise Exception("This is a test of exception handling.")
