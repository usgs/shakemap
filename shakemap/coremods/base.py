#stdlib imports
import logging

#local imports
from shakemap.utils.config import get_logging_config

class CoreModule(object):
    """Base class for any module in coremods which gets called by the shake program.

    

    """
    command_name = ''
    def __init__(self,eventid,testing=False):
        """Instantiate a CoreModule class with an event ID.

        """
        self._eventid = eventid
        self._testing = testing
        log_config = get_logging_config()
        log_name = log_config['loggers'].keys()[0]
        self.logger = logging.getLogger(log_name)

    def execute(self):
        pass

    def getDocString(self):
        return ''
