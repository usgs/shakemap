#stdlib imports
import logging
from abc import ABC, abstractmethod

#local imports
from shakemap.utils.config import get_logging_config

class CoreModule(ABC):
    """
    Base class for any module in coremods which gets called by the shake program.
    """
    command_name = ''
    def __init__(self, eventid):
        """Instantiate a CoreModule class with an event ID.

        """
        self._eventid = eventid
        log_config = get_logging_config()
        log_name = log_config['loggers'].keys()[0]
        self.logger = logging.getLogger(log_name)

    @abstractmethod
    def execute(self):
        pass
