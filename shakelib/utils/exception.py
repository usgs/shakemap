class ShakeLibException(Exception):
    """
    Class to represent errors in the Fault class.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
