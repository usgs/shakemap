# stdlib imports
import os.path
import json
from datetime import datetime

# third party imports
from shakelib.utils.containers import ShakeMapOutputContainer
from configobj import ConfigObj
from impactutils.transfer.factory import get_sender_class

# local imports
from .base import CoreModule
from shakemap.utils.config import get_config_paths,get_data_path

TIMEFMT = '%Y-%m-%d %H:%M:%S'

class TransferModule(CoreModule):
    """
    info -- Extract info.json from shake_result.hdf and write it as a file.
    """

    command_name = 'transfer'

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
            transfer_conf = os.path.join(install_path,'config','transfer.conf')
            if not os.path.isfile(transfer_conf):
                raise FileNotFoundError('%s does not exist.' % transfer_conf)

        # get the config information for transfer
        config = ConfigObj(transfer_conf,configspec=configspec)

        # get the output container with all the things in it
        datafile = os.path.join(products_dir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the ShakeMapOutputContainer and extract the data
        container = ShakeMapOutputContainer.load(datafile)

        # call the transfer method
        self._transfer(config, container, products_dir)

    def _transfer(self,config,container,products_dir):
        # extract the info.json object from the container
        info = json.loads(container.getString('info.json'))
        properties,product_properties = _get_properties(info)

        # get the config information:
        for transfer_method,transfer_dict in config.items():
            self.logger.debug('Doing %s transfers...' % transfer_method)
            for destination,params in transfer_dict.items():
                params.update(properties)
                fmt = 'Doing %s transfer to %s...'
                tpl = (transfer_method,destination)
                self.logger.debug(fmt % tpl)
                # append eventid to remote_directory, if present
                if 'remote_directory' in params:
                    rem_dir = params['remote_directory']
                    params['remote_directory'] = os.path.join(rem_dir,params['code'])

                # get the sender class from the type
                sender_class = get_sender_class(transfer_method)

                # try to send the data
                try:
                    if transfer_method == 'pdl':
                        sender = sender_class(properties = params,
                                              local_directory = products_dir,
                                              product_properties = product_properties)
                    else:
                        sender = sender_class(properties=params,
                                              local_directory = products_dir)
                    sender.send()
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
    product_properties = {}
    # origin info
    origin = info['input']['event_information']
    properties['eventsource'] = origin['eventsource']
    properties['eventsourcecode'] = origin['eventsourcecode']
    properties['code'] = origin['productcode']
    properties['source'] = origin['productsource']
    properties['type'] = origin['producttype']
    
    properties['magnitude'] = float(origin['magnitude'])
    properties['latitude'] = float(origin['latitude'])
    properties['longitude'] = float(origin['longitude'])
    properties['depth'] = float(origin['depth'])
    properties['eventtime'] = datetime.strptime(origin['origin_time'],TIMEFMT)

    
    product_properties['event-type'] = origin['event_type']
    product_properties['event-description'] = origin['event_description']

    # other metadata
    mmi_info = info['output']['ground_motions']['MMI']
    product_properties['maxmmi'] = mmi_info['max']
    product_properties['maxmmi-grid'] = mmi_info['max_grid']

    pgv_info = info['output']['ground_motions']['PGV']
    product_properties['maxpgv'] = pgv_info['max']
    product_properties['maxpgv-grid'] = pgv_info['max_grid']

    pga_info = info['output']['ground_motions']['PGA']
    product_properties['maxpga'] = pga_info['max']
    product_properties['maxpga-grid'] = pga_info['max_grid']

    psa03_info = info['output']['ground_motions']['SA(0.3)']
    product_properties['maxpsa03'] = psa03_info['max']
    product_properties['maxpsa03-grid'] = psa03_info['max_grid']

    psa10_info = info['output']['ground_motions']['SA(1.0)']
    product_properties['maxpsa10'] = psa10_info['max']
    product_properties['maxpsa10-grid'] = psa10_info['max_grid']

    psa30_info = info['output']['ground_motions']['SA(3.0)']
    product_properties['maxpsa30'] = psa30_info['max']
    product_properties['maxpsa30-grid'] = psa30_info['max_grid']

    product_properties['minimum-longitude'] = info['output']['map_information']['min']['longitude']
    product_properties['minimum-latitude'] = info['output']['map_information']['min']['latitude']
    product_properties['maximum-longitude'] = info['output']['map_information']['max']['longitude']
    product_properties['maximum-latitude'] = info['output']['map_information']['max']['latitude']

    product_properties['process-timestamp'] = info['processing']['shakemap_versions']['process_time']
    product_properties['version'] = info['processing']['shakemap_versions']['map_version']
    product_properties['map-status'] = info['processing']['shakemap_versions']['map_status']

    return (properties,product_properties)
                
        
