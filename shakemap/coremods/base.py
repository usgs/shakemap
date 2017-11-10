
class CoreModule(object):
    """Base class for any module in coremods which gets called by the shake program.

    

    """
    command_name = ''
    def __init__(self,eventid):
        """Instantiate a CoreModule class with an event ID.

        """
        self._eventid = eventid

    def execute(self):
        pass

    def getDocString(self):
        return ''
