#!/usr/bin/env python

# stdlib imports
import os.path

# local imports
from shakemap.utils.exception import ShakeMapException


class Sender(object):
    '''Base class for concrete subclasses that wrap around different methods of transmitting files.
    '''

    def __init__(self, properties=None, files=None, directory=None):
        self.properties = properties
        if files is not None:
            if not isinstance(files, list):
                raise ShakeMapException('Input files must be a list')
            for f in files:
                if not os.path.isfile(f):
                    raise ShakeMapException(
                        'Input file %s could not be found' % f)
        if directory is not None:
            if not os.path.isdir(directory):
                raise ShakeMapException(
                    'Input directory %s could not be found' % directory)
        if files is not None:
            self.files = files
        else:
            self.files = []
        if directory is not None:
            self.directories = [directory]
        else:
            self.directories = []

    def addProperty(self, key, value):
        self.properties[key] = value

    def addFiles(self, files):
        for f in files:
            if not os.path.isfile(f):
                raise ShakeMapException('Input file %s could not be found' % f)
        self.files += files

    def addDirectory(self, directory):
        if not os.path.isdir(directory):
            raise ShakeMapException(
                'Input directory %s could not be found' % directory)
        self.directories += directory

    # this is implemented in the subclasses
    def send(self):
        pass

    # this is implemented in the subclasses
    def delete(self):
        pass
