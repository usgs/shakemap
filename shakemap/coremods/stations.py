#stdlib imports
import sys
import os.path
import json
import logging

#third party imports
from shakelib.utils.containers import OutputContainer
from configobj import ConfigObj

#local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths,get_logging_config

ALLLOWED_FORMATS = ['json']

class StationModule(CoreModule):
    """
    **stations** -- Generate stationlist.json from shake_result.hdf.
    """
    command_name = 'stations'
    def execute(self):
        """Write stationlist.json file.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not exist.
        """
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        datafile = os.path.join(datadir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)
            
        # Open the OutputContainer and extract the data
        container = OutputContainer.load(datafile)

        # get the path to the products.conf file, load the config
        config_file = os.path.join(install_path, 'config', 'products.conf')
        config = ConfigObj(config_file)

        # create ShakeMap station data file
        formats = config['products']['info']['formats']
        for fformat in formats:
            if fformat not in ALLLOWED_FORMATS:
                logger.warn('Specified format %s not in list of defined formats.  Skipping.' % fformat)
                continue
            if fformat == 'json':
                self.logger.info('Writing rupture.json file...')
                stationstring = container.getString('stationlist.json')
                station_file = os.path.join(datadir,'stationlist.json')
                f = open(station_file,'w')
                f.write(stationstring)
                f.close()
        

