# stdlib imports
import os.path
import json

# third party imports
from shakelib.utils.containers import ShakeMapOutputContainer
from configobj import ConfigObj
from impactutils.transfer.factory import get_sender_class

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths

class TransferModule(CoreModule):
    """
    info -- Extract info.json from shake_result.hdf and write it as a file.
    """

    command_name = 'info'

    def execute(self):
        """
        Tranfer ShakeMap products using methods configured in transfer.conf.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        products_dir = os.path.join(datadir, 'products')
        if not os.path.isdir(products_dir):
            raise NotADirectoryError('%s does not exist.' % products_dir)

        # get the path to the transfer.conf spec file
        configspec = os.path.join(get_data_path(),'transferspec.conf')
        
        # look for an event specific transfer.conf file
        transfer_conf = os.path.join(datadir,'transfer.conf')
        if not os.path.isfile(transfer_conf):
            # if not there, use the system one
            transfer_conf = os.path.join(install_path,'transfer.conf')
            if not os.path.isfile(transfer_conf):
                raise FileNotFoundError('%s does not exist.' % transfer_conf)

        # get the config information for transfer
        config = ConfigObj(transfer_conf,configspec=configspec)

        # get the output container with all the things in it
        datafile = os.path.join(datadir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the ShakeMapOutputContainer and extract the data
        container = ShakeMapOutputContainer.load(datafile)

        # call the transfer method
        transfer(config, container, products_dir)

def _transfer(config,container,products_dir):
    # extract the info.json object from the container
    info = json.loads(container.getString('info.json'))
    properties = _get_properties(info)

    # get the config information:
    for transfer_method,transfer_dict in config.items():
        self.logger.debug('Doing %s transfers...' % transfer_method)
        for destination,params in transfer_dict.items():
            fmt = 'Doing %s transfer to %s...'
            tpl = (transfer_method,destination)
            self.logger.debug(fmt % tpl)
            # append eventid to remote_directory, if present
            if 'remote_directory' in params:
                rem_dir = params['remote_directory']
                params['remote_directory'] = os.path.join(rem_dir,self._eventid)
            params['code'] = self._eventid
            params['eventsource'] = eventsource
            params['eventsourcecode'] = self._eventid.replace(eventsource,'')
            params['magnitude'] = magnitude
            params['latitude'] = latitude
            params['longitude'] = longitude
            params['depth'] = depth
            params['eventtime'] = eventtime

            # get the sender class from the type
            sender_class = get_sender_class(transfer_method)

            # try to send the data
            try:
                if method == 'pdl':
                    sender = sender_class(properties = params,
                                          local_directory = products_dir,
                                          product_properties = properties)
                else:
                    sender = sender_class(properties=params,
                                          local_directory = products_dir)

            except Exception as e:
                # for the standard config, this should generate an email to the developer
                # list.
                msg = str(e)
                fmt = 'Transfer for %s method, %s destination failed with error "%s".'
                tpl = (transfer_method,destination,msg)
                self.logger.exception(fmt % tpl)
                continue
                
def _get_properties(info):
    properties = {}
    # origin info
    origin = info['input']['event_information']
    properties['eventsource'] = origin['eventsource']
    properties['magnitude'] = origin['magnitude']
    properties['latitude'] = origin['latitude']
    properties['longitude'] = origin['longitude']
    properties['depth'] = origin['depth']
    properties['eventtime'] = origin['origin_time']
    properties['event-type'] = origin['event_type']
    properties['event-description'] = origin['event_description']

    # other metadata
    mmi_info = info['output']['ground_motions']['mmi']
    properties['maxmmi'] = mmi_info['max']
    properties['maxmmi-grid'] = mmi_info['max_grid']

    pgv_info = info['output']['ground_motions']['pgv']
    properties['maxpgv'] = pgv_info['max']
    properties['maxpgv-grid'] = pgv_info['max_grid']

    pga_info = info['output']['ground_motions']['pga']
    properties['maxpga'] = pga_info['max']
    properties['maxpga-grid'] = pga_info['max_grid']

    psa03_info = info['output']['ground_motions']['psa03']
    properties['maxpsa03'] = psa03_info['max']
    properties['maxpsa03-grid'] = psa03_info['max_grid']

    psa10_info = info['output']['ground_motions']['psa10']
    properties['maxpsa10'] = psa10_info['max']
    properties['maxpsa10-grid'] = psa10_info['max_grid']

    psa30_info = info['output']['ground_motions']['psa30']
    properties['maxpsa30'] = psa30_info['max']
    properties['maxpsa30-grid'] = psa30_info['max_grid']

    properties['minimum-longitude'] = info['output']['map_information']['min']['longitude']
    properties['minimum-latitude'] = info['output']['map_information']['min']['latitude']
    properties['maximum-longitude'] = info['output']['map_information']['max']['longitude']
    properties['maximum-latitude'] = info['output']['map_information']['max']['latitude']

    properties['process-timestamp'] = info['processing']['shakemap_versions']['process_time']
    properties['version'] = info['processing']['shakemap_versions']['map_version']
    properties['map-status'] = info['processing']['shakemap_versions']['map_status']

    return properties
                
        
