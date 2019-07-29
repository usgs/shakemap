class TerminateShakeMap(Exception):
    """
    Class to terminate ShakeMap processing without an error.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
