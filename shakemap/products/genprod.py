#!/usr/bin/env python

#stdlib imports
import json
from datetime import datetime

from mapio.shake import ShakeGrid
import shakemap

TIMEFMT = '%Y-%m-%d %H:%M:%S'

def make_xml_grid(container,xml_type='grid'):
    """
    Create an XML file containing either the ground motion arrays from sm_model,
    or the uncertainty values for those same ground motions.
    
    Args:
        container (GridHDFContainer): The container created by sm_model.
        xml_type (str): A string, either 'grid' or 'uncertainty'.
    Returns:
        (ShakeGrid): A ShakeGrid data structure containing ground motion or uncertain
        arrays, and associated metadata.
    """

    #get all of the grid layers and the geodict
    gridnames = container.getGridNames()
    layers = {}
    field_keys = {}
    for gridname in gridnames:
        if xml_type == 'grid' and gridname.endswith('_sd'):
            continue
        if xml_type == 'uncertainty' and not gridname.endswith('_sd'):
            continue
        grid,metadata = container.getGrid(gridname)
        layers[gridname] = grid.getData()
        units = metadata['units']
        digits = metadata['digits']
        field_keys[gridname] = (units,digits)
    geodict = grid.getGeoDict()

    config = container.getMetadata('config')

    #event dictionary
    info_data,info_metadata = container.getData('info.json')
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

    return shake_grid
    

