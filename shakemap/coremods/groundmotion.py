# stdlib imports
import os.path
from io import StringIO
import json

# third party imports
from libcomcat.search import get_event_by_id
from libcomcat.classes import DetailEvent
from impactutils.io.table import dataframe_to_xml
import pandas as pd
import numpy as np

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths

# Get rid of stupid pandas warning
pd.options.mode.chained_assignment = None

HOST = 'dev-earthquake.cr.usgs.gov'


class GroundMotionModule(CoreModule):
    """
    groundmotion -- Search ComCat for GroundMotion data and turn it into a ShakeMap data file.
    """

    command_name = 'groundmotion'

    def execute(self):
        """
        Write ugroundmotion_dat.json metadata file.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        _, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        if not os.path.isdir(datadir):
            os.makedirs(datadir)

        # try to find the event by our event id
        try:
            detail = get_event_by_id(self._eventid, host=HOST)
            if not detail.hasProduct('ground-motion'):
                return
            groundmotion = detail.getProducts('ground-motion')[0]
            fname = 'groundmotions_dat.json'
            gbytes, gurl = groundmotion.getContentBytes(fname)
            outname = os.path.join(datadir, 'ugroundmotions_dat.json')
            with open(outname, 'wt') as f:
                f.write(gbytes.decode('utf-8'))
            self.logger.info('Created ground motions data file %s' % outname)
        except Exception as e:
            fmt = 'Could not retrieve ground motion data for %s - error "%s"'
            self.logger.warning(fmt % (self._eventid, str(e)))
            return
