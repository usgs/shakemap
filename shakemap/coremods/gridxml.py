#stdlib imports
import sys
import os.path
import json
import logging
import glob
import warnings
from datetime import datetime

#third party imports
from shapely.geometry import MultiLineString,mapping
import fiona
import numpy as np
from skimage import measure
from scipy.ndimage.filters import median_filter
from shakelib.utils.containers import OutputContainer
from mapio.shake import ShakeGrid

#local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths
import shakemap

# historically, we only had the component we are now calling 'Larger'.
# As we do not intend grid.xml files to be forward compatible with
# additional layers of information and different components (rotd50, etc.)
# we'll hard code this here until grid.xml files experience their heat death.
COMPONENT = 'Larger'

TIMEFMT = '%Y-%m-%d %H:%M:%S'

class GridXMLModule(CoreModule):
    """contour - Generate contours of all configured IMT values.
    """
    command_name = 'gridxml'
    def execute(self):
        logger = logging.getLogger(__name__)
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            print('%s is not a valid directory.' % datadir)
            sys.exit(1)
        datafile = os.path.join(datadir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            print('%s is not a valid shake result file.' % datafile)
            sys.exit(1)
        
        # Open the OutputContainer and extract the data
        container = OutputContainer.load(datafile)

        #get all of the grid layers and the geodict
        gridnames = container.getIMTs(COMPONENT)
        layers = {}
        field_keys = {}
        xml_types = ['grid','uncertainty']
        for xml_type in xml_types:
            for gridname in gridnames:
                imtdict = container.getIMT(gridname,COMPONENT)
                if xml_type == 'grid':
                    grid = imtdict['mean']
                    metadata = imtdict['mean_metadata']
                elif xml_type == 'uncertainty':
                    grid = imtdict['mean']
                    metadata = imtdict['mean_metadata']

                units = metadata['units']
                digits = metadata['digits']
                grid_data = grid.getData()
                #convert from HDF units to legacy grid.xml units
                if units == 'ln(cm/s)':
                    grid_data = np.exp(grid_data)
                    units = 'cm/s'
                elif units == 'ln(g)':
                    grid_data = np.exp(grid_data)*100
                    units = '%g'
                else:
                    pass
                layers[gridname] = grid_data

                field_keys[gridname] = (units,digits)
            geodict = grid.getGeoDict()

            config = container.getDictionary('config')

            #event dictionary
            info_data = container.getString('info.json')
            info = json.loads(info_data)
            event_info = info['input']['event_information']
            event_dict = {}
            event_dict['event_id'] = event_info['event_id']
            event_dict['magnitude'] = float(event_info['magnitude'])
            event_dict['depth'] = float(event_info['depth'])
            event_dict['lat'] = float(event_info['latitude'])
            event_dict['lon'] = float(event_info['longitude'])
            event_dict['event_timestamp'] = datetime.strptime(event_info['origin_time'],TIMEFMT)
            event_dict['event_description'] = event_info['location']
            #TODO the following is SUPER-SKETCHY - we need to save the event network info!!!
            event_dict['event_network'] = event_dict['event_id'][0:2]

            #shake dictionary
            shake_dict = {}
            shake_dict['event_id'] = event_dict['event_id']
            shake_dict['shakemap_id'] = event_dict['event_id']
            #TODO - where are we supposed to get shakemap version
            shake_dict['shakemap_version'] = 1
            shake_dict['code_version'] = shakemap.__version__
            shake_dict['process_timestamp'] = datetime.utcnow()
            shake_dict['shakemap_originator'] = config['system']['source_network']
            shake_dict['map_status'] = config['system']['map_status']
            #TODO - we need a source for this!!!
            shake_dict['shakemap_event_type'] = 'ACTUAL'

            shake_grid = ShakeGrid(layers,geodict,event_dict,shake_dict,{},field_keys=field_keys)
            fname = os.path.join(datadir,'%s.xml' % xml_type)
            logger.info('Saving IMT grids to %s' % fname)
            shake_grid.save(fname) #TODO - set grid version number
