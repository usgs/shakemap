"""
Reads shake_data.hdf from the event's current directory and adds local
configs, data, etc., then writes a new shake_data.hdf.
"""

# stdlib imports
import argparse
import inspect
import os.path
import glob
import datetime
import shutil
import sys
import logging

# third party imports
from configobj import ConfigObj

# local imports
from .base import CoreModule
from shakelib.utils.containers import ShakeMapInputContainer
from shakemap.utils.config import (get_config_paths,
                                   get_custom_validator,
                                   config_error,
                                   check_config,
                                   get_configspec,
                                   path_macro_sub)
from shakemap.utils.logging import get_logging_config


class AugmentModule(CoreModule):
    """
    augment -- Incorporate additional content into the shake_data.hdf input
                     file.
    """

    command_name = 'augment'

    def __init__(self, eventid, comment=None):
        """
        Instantiate a CoreModule class with an event ID.
        """
        super(AugmentModule, self).__init__(eventid)
        if comment is not None:
            self.comment = comment

    def execute(self):
        """
        Augment a ShakeMap input data file with local configs, data, rupture,
        etc. The version history will only be incremented if the originator
        differs from the originator in the previous line of the history.

        Raises:
            NotADirectoryError: When the event data directory does not
                exist.
            FileNotFoundError: When the the event's event.xml file does
                not exist.
            RuntimeError: When there are problems parsing the configuration.
            ValidateError: When there are configuration items missing or mis-
                configured.
        """

        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)

        hdf_file = os.path.join(datadir, 'shake_data.hdf')
        if not os.path.isfile(hdf_file):
            raise FileNotFoundError('%s does not exist. Use assemble.' %
                                    hdf_file)
        shake_data = ShakeMapInputContainer.load(hdf_file)

        # Prompt for a comment string if none is provided on the command line
        if self.comment is None:
            if sys.stdout is not None and sys.stdout.isatty():
                self.comment = input(
                        'Please enter a comment for this version.\n'
                        '(Start with "+" if you wish to append to the\n'
                        'existing comment; "+" by itself will preserve\n'
                        'existing comments.)\n'
                        'comment: ')
            else:
                self.comment = ''

        #
        # Clear away results from previous runs
        #
        products_path = os.path.join(datadir, 'products')
        if os.path.isdir(products_path):
            shutil.rmtree(products_path, ignore_errors=True)
        pdl_path = os.path.join(datadir, 'pdl')
        if os.path.isdir(pdl_path):
            shutil.rmtree(pdl_path, ignore_errors=True)

        #
        # Get the config from the HDF file and merge in the local configs
        #
        spec_file = get_configspec()
        validator = get_custom_validator()
        shake_config = shake_data.getConfig()
        shake_config = ConfigObj(shake_config, configspec=spec_file)
        #
        # This is a weird hack to get around a bug/feature of ConfigObj
        # that results in the validation failing if max_workers is already
        # an integer.
        #
        if 'max_workers' in shake_config['system']:
            shake_config['system']['max_workers'] = \
                    str(shake_config['system']['max_workers'])

        modules_file = os.path.join(install_path, 'config', 'modules.conf')
        if os.path.isfile(modules_file):
            self.logger.debug('Found a modules file.')
            modules = ConfigObj(modules_file, configspec=spec_file)
            shake_config.merge(modules)
        gmpe_file = os.path.join(install_path, 'config', 'gmpe_sets.conf')
        if os.path.isfile(gmpe_file):
            self.logger.debug('Found a gmpe file.')
            gmpe_sets = ConfigObj(gmpe_file, configspec=spec_file)
            shake_config.merge(gmpe_sets)
        config_file = os.path.join(install_path, 'config', 'model.conf')
        if os.path.isfile(config_file):
            self.logger.debug('Found a global config file.')
            global_config = ConfigObj(config_file, configspec=spec_file)
            shake_config.merge(global_config)

        # extent conf (may not be present)
        extent_config = os.path.join(install_path, 'config', 'extent.conf')
        if os.path.isfile(extent_config):
            extent_config = ConfigObj(extent_config,
                                      configspec=spec_file)
        else:
            extent_config = ConfigObj()
        shake_config.merge(extent_config)
        #
        # this is the event specific model.conf (may not be present)
        # prefer model.conf to model_select.conf
        #
        event_config_file = os.path.join(datadir, 'model.conf')
        event_config_zc_file = os.path.join(datadir, 'model_select.conf')
        if os.path.isfile(event_config_file):
            self.logger.debug('Found an event specific model.conf file.')
            event_config = ConfigObj(event_config_file,
                                     configspec=spec_file)
            shake_config.merge(event_config)
        elif os.path.isfile(event_config_zc_file):
            self.logger.debug('Found an event specific model_select file.')
            event_config = ConfigObj(event_config_zc_file,
                                     configspec=spec_file)
            shake_config.merge(event_config)
        #
        # Validate the resulting config
        #
        results = shake_config.validate(validator)
        if not results or isinstance(results, dict):
            config_error(shake_config, results)
        check_config(shake_config, self.logger)

        global_data_path = os.path.join(os.path.expanduser('~'),
                                        'shakemap_data')
        #
        # If there is a prediction_location->file file, then we need
        # to expand any macros
        #
        if 'file' in shake_config['interp']['prediction_location']:
            loc_file = shake_config['interp']['prediction_location']['file']
            if loc_file and loc_file != 'None':      # 'None' is a string here
                loc_file = path_macro_sub(loc_file, ip=install_path,
                                          dp=data_path, gp=global_data_path,
                                          ei=self._eventid)
                if not os.path.isfile(loc_file):
                    raise FileNotFoundError('prediction file "%s" is not a '
                                            'valid file' % loc_file)
                shake_config['interp']['prediction_location']['file'] = \
                    loc_file
        #
        # Put the updated config back into shake_data.hdf`
        #
        config = shake_config.dict()
        shake_data.setConfig(config)
        #
        # Look for additional data files and update the stationlist if found
        #
        datafiles = glob.glob(os.path.join(datadir, '*_dat.xml'))
        if os.path.isfile(os.path.join(datadir, 'stationlist.xml')):
            datafiles.append(os.path.join(datadir, 'stationlist.xml'))
        datafiles += glob.glob(os.path.join(datadir, '*_dat.json'))
        if os.path.isfile(os.path.join(datadir, 'stationlist.json')):
            datafiles.append(os.path.join(datadir, 'stationlist.json'))
        if datafiles:
            self.logger.debug('Found additional data files...')
            shake_data.addStationData(datafiles)
        #
        # Look for a rupture file and replace the existing one if found
        #
        rupturefile = os.path.join(datadir, 'rupture.json')
        eventxml = os.path.join(datadir, 'event.xml')
        if not os.path.isfile(eventxml):
            eventxml = None
        if not os.path.isfile(rupturefile):
            faultfiles = glob.glob(os.path.join(datadir, '*_fault.txt'))
            if len(faultfiles):
                rupturefile = faultfiles[0]
            else:
                rupturefile = None
        if (rupturefile and os.path.isfile(rupturefile)) \
                or eventxml is not None:
            self.logger.debug('Updating rupture/origin information.')
            shake_data.updateRupture(
                eventxml=eventxml, rupturefile=rupturefile)

        #
        # Sort out the version history. We're working with an existing
        # HDF file, so: if we are the originator, just update the timestamp,
        # otherwise add a new line.
        #
        timestamp = datetime.datetime.utcnow().strftime('%FT%TZ')
        originator = config['system']['source_network']

        history = shake_data.getVersionHistory()
        if history['history'][-1][1] == originator:
            history['history'][-1][0] = timestamp
            if self.comment.startswith('+'):
                if self.comment.replace('+', '') != '':
                    history['history'][-1][3] += self.comment.replace('+', ' ')
            else:
                history['history'][-1][3] = self.comment
        else:
            version = int(history['history'][-1][2]) + 1
            if self.comment.startswith('+'):
                new_line = [timestamp, originator, version,
                            self.comment.replace('+', '')]
            else:
                new_line = [timestamp, originator, version, self.comment]
            history['history'].append(new_line)
        shake_data.setVersionHistory(history)

        shake_data.close()

    def parseArgs(self, arglist):
        """
        Set up the object to accept the --comment flag.
        """
        parser = argparse.ArgumentParser(
            prog=self.__class__.command_name,
            description=inspect.getdoc(self.__class__))
        parser.add_argument('-c', '--comment', help='Provide a comment for '
                            'this version of the ShakeMap. If the comment '
                            'has spaces, the string should be quoted (e.g., '
                            '--comment "This is a comment.")')
        #
        # This line should be in any modules that overrides this
        # one. It will collect up everything after the current
        # modules options in args.rem, which should be returned
        # by this function. Note: doing parser.parse_known_args()
        # will not work as it will suck up any later modules'
        # options that are the same as this one's.
        #
        parser.add_argument('rem', nargs=argparse.REMAINDER,
                            help=argparse.SUPPRESS)
        args = parser.parse_args(arglist)
        self.comment = args.comment
        return args.rem
