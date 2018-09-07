# stdlib imports
import os.path
from datetime import datetime
import glob
import re
import shutil
import logging

# third party imports
from impactutils.io.smcontainers import ShakeMapOutputContainer
from configobj import ConfigObj
from validate import Validator
from impactutils.transfer.factory import get_sender_class

# local imports
from .base import CoreModule
from shakemap.utils.config import (get_config_paths,
                                   get_data_path,
                                   config_error)
from shakemap.utils.macros import get_macros
from shakelib.rupture import constants

NO_TRANSFER = 'NO_TRANSFER'

# what are the names of the cancel files for the different transfer methods
# that use them?
CANCEL_FILES = {'ftp': 'CANCEL',
                'ssh': 'CANCEL',
                'email': None,
                'copy': 'CANCEL'}


class TransferModule(CoreModule):
    """
    transfer -- Transfer data via any configured methods (ftp, pdl, etc.)
    """

    command_name = 'transfer'
    dependencies = [('products/*', False)]

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

        # look for the presence of a NO_TRANSFER file in the datadir.
        notransfer = os.path.join(datadir, NO_TRANSFER)
        if os.path.isfile(notransfer):
            self.logger.info(
                'Event has a %s file blocking transfer.' % NO_TRANSFER)
            return

        pdl_dir = os.path.join(datadir, 'pdl')
        products_dir = os.path.join(datadir, 'products')
        if not os.path.isdir(pdl_dir):
            raise NotADirectoryError('%s does not exist.' % pdl_dir)

        # get the path to the transfer.conf spec file
        configspec = os.path.join(get_data_path(), 'transferspec.conf')

        # look for an event specific transfer.conf file
        transfer_conf = os.path.join(datadir, 'transfer.conf')
        if not os.path.isfile(transfer_conf):
            # if not there, use the system one
            transfer_conf = os.path.join(
                install_path, 'config', 'transfer.conf')
            if not os.path.isfile(transfer_conf):
                raise FileNotFoundError('%s does not exist.' % transfer_conf)

        # get the config information for transfer
        config = ConfigObj(transfer_conf, configspec=configspec)
        results = config.validate(Validator())
        if not isinstance(results, bool) or not results:
            config_error(config, results)

        # get the output container with all the things in it
        datafile = os.path.join(products_dir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the ShakeMapOutputContainer and extract the data
        container = ShakeMapOutputContainer.load(datafile)
        # extract the info.json object from the container
        info = container.getMetadata()
        container.close()

        # call the transfer method
        _transfer(config, info, pdl_dir, products_dir)

        # copy the current folder to a new backup directory
        self._make_backup(data_path)

    def _make_backup(self, data_path):
        data_dir = os.path.join(data_path, self._eventid)
        current_dir = os.path.join(data_dir, 'current')
        backup_dirs = glob.glob(os.path.join(data_dir, 'backup*'))
        latest_version = 0

        # and get the most recent version number
        for backup_dir in backup_dirs:
            if not os.path.isdir(backup_dir):
                continue
            match = re.search('[0-9]*$', backup_dir)
            if match is not None:
                version = int(match.group())
                if version > latest_version:
                    latest_version = version

        new_version = latest_version + 1
        backup = os.path.join(data_dir, 'backup%04i' % new_version)
        shutil.copytree(current_dir, backup)
        self.logger.debug('Created backup directory %s' % backup)


def _transfer(config, info, pdl_dir, products_dir, cancel=False):

    # get the macros that may be in the email sender config
    macros = get_macros(info)

    properties, product_properties = _get_properties(info)

    # get the config information:
    for transfer_method, transfer_dict in config.items():
        logging.debug('Doing %s transfers...' % transfer_method)
        for destination, params in transfer_dict.items():
            params.update(properties)
            fmt = 'Doing %s transfer to %s...'
            tpl = (transfer_method, destination)
            logging.debug(fmt % tpl)
            # append eventid to remote_directory, if present
            if 'remote_directory' in params:
                rem_dir = params['remote_directory']
                params['remote_directory'] = os.path.join(
                    rem_dir, params['code'])

            # get the sender class from the type
            sender_class = get_sender_class(transfer_method)

            # try to send the data
            try:
                if transfer_method == 'pdl':
                    downloads_dir = os.path.join(pdl_dir, 'download')
                    if os.path.isdir(downloads_dir):
                        shutil.rmtree(downloads_dir, ignore_errors=True)
                    shutil.copytree(products_dir, downloads_dir)
                    sender = sender_class(
                        properties=params,
                        local_directory=pdl_dir,
                        product_properties=product_properties)
                elif transfer_method == 'email':
                    # replace macro strings with actual strings
                    for pkey, param in params.items():
                        for macro, replacement in macros.items():
                            if isinstance(param, str):
                                param = param.replace('[%s]' % macro,
                                                      replacement)
                                params[pkey] = param

                    # get full path to all file attachments
                    attachments = []
                    for lfile in params['attachments']:
                        fullfile = os.path.join(products_dir, lfile)
                        if not os.path.isfile(fullfile):
                            logging.warn('%s does not exist.' % fullfile)
                        attachments.append(fullfile)

                    # construct an instance of the EmailSender
                    sender = sender_class(properties=params,
                                          local_files=attachments)
                else:
                    cancelfile = CANCEL_FILES[transfer_method]
                    sender = sender_class(properties=params,
                                          local_directory=products_dir,
                                          cancelfile=cancelfile)

                if cancel:
                    msg = sender.cancel()
                else:
                    nfiles, msg = sender.send()
                fmt = '%i files sent.  Message from sender: \n"%s"'
                tpl = (nfiles, msg)
                logging.info(fmt % tpl)
                if transfer_method == 'pdl':
                    shutil.rmtree(downloads_dir, ignore_errors=True)
            except Exception as e:
                # for the standard config, this should generate an email to
                # the developer list.
                msg = str(e)
                fmt = ('Transfer for %s method, %s destination '
                       'failed with error "%s".')
                tpl = (transfer_method, destination, msg)
                logging.warning(fmt % tpl)
                continue


def _get_properties(info):
    properties = {}
    product_properties = {}
    # origin info
    origin = info['input']['event_information']
    properties['eventsource'] = origin['netid']
    # The netid could be a valid part of the eventsourcecode, so we have to
    # check here if it ***starts with*** the netid
    # This fix should already be done by the time we get here, but this
    # is just an insurance check
    if origin['eventsourcecode'].startswith(origin['netid']):
        eventsourcecode = origin['eventsourcecode'].replace(origin['netid'],
                                                            '', 1)
    else:
        eventsourcecode = origin['eventsourcecode']
    properties['eventsourcecode'] = eventsourcecode
    properties['code'] = origin['productcode']
    properties['source'] = origin['productsource']
    properties['type'] = origin['producttype']

    properties['magnitude'] = float(origin['magnitude'])
    properties['latitude'] = float(origin['latitude'])
    properties['longitude'] = float(origin['longitude'])
    properties['depth'] = float(origin['depth'])
    try:
        properties['eventtime'] = datetime.strptime(origin['origin_time'],
                                                    constants.TIMEFMT)
    except ValueError:
        properties['eventtime'] = datetime.strptime(origin['origin_time'],
                                                    constants.ALT_TIMEFMT)

    product_properties['event-type'] = origin['event_type']
    product_properties['event-description'] = origin['event_description']

    # other metadata
    if 'MMI' in info['output']['ground_motions']:
        mmi_info = info['output']['ground_motions']['MMI']
        product_properties['maxmmi'] = mmi_info['max']
        product_properties['maxmmi-grid'] = mmi_info['max_grid']

    if 'PGV' in info['output']['ground_motions']:
        pgv_info = info['output']['ground_motions']['PGV']
        product_properties['maxpgv'] = pgv_info['max']
        product_properties['maxpgv-grid'] = pgv_info['max_grid']

    if 'PGA' in info['output']['ground_motions']:
        pga_info = info['output']['ground_motions']['PGA']
        product_properties['maxpga'] = pga_info['max']
        product_properties['maxpga-grid'] = pga_info['max_grid']

    if 'SA(0.3)' in info['output']['ground_motions']:
        psa03_info = info['output']['ground_motions']['SA(0.3)']
        product_properties['maxpsa03'] = psa03_info['max']
        product_properties['maxpsa03-grid'] = psa03_info['max_grid']

    if 'SA(1.0)' in info['output']['ground_motions']:
        psa10_info = info['output']['ground_motions']['SA(1.0)']
        product_properties['maxpsa10'] = psa10_info['max']
        product_properties['maxpsa10-grid'] = psa10_info['max_grid']

    if 'SA(3.0)' in info['output']['ground_motions']:
        psa30_info = info['output']['ground_motions']['SA(3.0)']
        product_properties['maxpsa30'] = psa30_info['max']
        product_properties['maxpsa30-grid'] = psa30_info['max_grid']

    product_properties['minimum-longitude'] = \
        info['output']['map_information']['min']['longitude']
    product_properties['minimum-latitude'] = \
        info['output']['map_information']['min']['latitude']
    product_properties['maximum-longitude'] = \
        info['output']['map_information']['max']['longitude']
    product_properties['maximum-latitude'] = \
        info['output']['map_information']['max']['latitude']

    product_properties['process-timestamp'] = \
        info['processing']['shakemap_versions']['process_time']
    product_properties['version'] = \
        info['processing']['shakemap_versions']['map_version']
    product_properties['map-status'] = \
        info['processing']['shakemap_versions']['map_status']

    return (properties, product_properties)
