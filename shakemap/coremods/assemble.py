"""
Collect configuration, station data, finite fault data, etc., into
an InputContainer and write it out as shake_data.hdf.
"""

# stdlib imports
import os.path
import glob
import datetime
import shutil

# third party imports
from configobj import ConfigObj

# local imports
from .base import CoreModule
from shakelib.utils.containers import ShakeMapInputContainer
from shakemap.utils.config import get_config_paths, get_custom_validator,\
    config_error, check_config, get_configspec


class AssembleModule(CoreModule):
    """
    assemble -- Assemble ShakeMap input data into the shake_data.hdf input
                      file.
    """

    command_name = 'assemble'

    def execute(self):
        """
        Assemble ShakeMap input data and write and ShakeMapInputContainer named
        shake_data.hdf in the event's 'current' directory.

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

        eventxml = os.path.join(datadir, 'event.xml')
        self.logger.debug('Looking for event.xml file...')
        if not os.path.isfile(eventxml):
            raise FileNotFoundError('%s does not exist.' % eventxml)

        # find any source.txt or moment.xml files
        momentfile = os.path.join(datadir, 'moment.xml')
        sourcefile = os.path.join(datadir, 'source.txt')
        if not os.path.isfile(sourcefile):
            sourcefile = None
        if not os.path.isfile(momentfile):
            momentfile = None

        #
        # Clear away results from previous runs
        #
        products_path = os.path.join(datadir, 'products')
        if os.path.isdir(products_path):
            shutil.rmtree(products_path, ignore_errors=True)

        #
        # Look for global configs in install_path/config
        #
        spec_file = get_configspec()
        validator = get_custom_validator()
        self.logger.debug('Looking for configuration files...')
        modules = ConfigObj(
            os.path.join(install_path, 'config', 'modules.conf'),
            configspec=spec_file)
        gmpe_sets = ConfigObj(
            os.path.join(install_path, 'config', 'gmpe_sets.conf'),
            configspec=spec_file)
        global_config = ConfigObj(
            os.path.join(install_path, 'config', 'model.conf'),
            configspec=spec_file)

        #
        # this is the event specific model.conf (may not be present)
        # prefer model.conf to model_zc.conf
        #
        event_config_file = os.path.join(datadir, 'model.conf')
        event_config_zc_file = os.path.join(datadir, 'model_zc.conf')
        if os.path.isfile(event_config_file):
            event_config = ConfigObj(event_config_file,
                                     configspec=spec_file)
        elif os.path.isfile(event_config_zc_file):
            event_config = ConfigObj(event_config_zc_file,
                                     configspec=spec_file)
        else:
            event_config = ConfigObj()

        #
        # start merging event_config
        #
        global_config.merge(event_config)
        global_config.merge(modules)
        global_config.merge(gmpe_sets)

        results = global_config.validate(validator)
        if not isinstance(results, bool) or not results:
            config_error(global_config, results)

        check_config(global_config, self.logger)

        #
        # The vs30 file may have macros in it
        #
        vs30file = global_config['data']['vs30file']
        if vs30file:
            vs30file = vs30file.replace('<INSTALL_DIR>', install_path)
            vs30file = vs30file.replace('<DATA_DIR>', data_path)
            vs30file = vs30file.replace('<EVENT_ID>', self._eventid)
            if not os.path.isfile(vs30file):
                raise FileNotFoundError("vs30 file '%s' is not a valid file" %
                                        vs30file)
            global_config['data']['vs30file'] = vs30file
        #
        # If there is a prediction_location->file file, then we need
        # to expand any macros
        #
        if 'file' in global_config['interp']['prediction_location']:
            loc_file = global_config['interp']['prediction_location']['file']
            if loc_file and loc_file != 'None':      # 'None' is a string here
                loc_file = loc_file.replace('<INSTALL_DIR>', install_path)
                loc_file = loc_file.replace('<DATA_DIR>', data_path)
                loc_file = loc_file.replace('<EVENT_ID>', self._eventid)
                if not os.path.isfile(loc_file):
                    raise FileNotFoundError("prediction file '%s' is not "
                                            "a valid file" % loc_file)
                global_config['interp']['prediction_location']['file'] = \
                    loc_file

        config = global_config.dict()

        self.logger.debug('Looking for data files...')
        datafiles = glob.glob(os.path.join(datadir, '*_dat.xml'))
        if os.path.isfile(os.path.join(datadir, 'stationlist.xml')):
            datafiles.append(os.path.join(datadir, 'stationlist.xml'))

        self.logger.debug('Looking for rupture files...')
        rupturefiles = glob.glob(os.path.join(datadir, '*_fault.txt'))
        rupturefile = None
        if len(rupturefiles):
            rupturefile = rupturefiles[0]
        #
        # Sort out the version history. Get the most recent backup file and
        # extract the existing history. Then add a new line for this run.
        #
        timestamp = datetime.datetime.utcnow().strftime('%FT%TZ')
        originator = config['system']['source_network']
        backup_dirs = sorted(
            glob.glob(os.path.join(datadir, '..', '.backup*')),
            reverse=True)
        if len(backup_dirs):
            #
            # Backup files exist so find the latest one and extract its
            # history, then add a new line that increments the version
            #
            bu_file = os.path.join(backup_dirs[0], 'shake_data.hdf')
            bu_ic = ShakeMapInputContainer.load(bu_file)
            history = bu_ic.getVersionHistory()
            bu_ic.close()
            version = int(
                backup_dirs[0].replace(
                    os.path.join(datadir, '..', '.backup'), ''))
            version += 1
            new_line = [timestamp, originator, version]
            history['history'].append(new_line)
        elif os.path.isfile(os.path.join(datadir, 'shake_data.hdf')):
            #
            # No backups are available, but there is an existing shake_data
            # file. Extract its history and update the timestamp and
            # source network (but leave the version alone).
            # If there is no history, just start a new one with version 1
            #
            bu_file = os.path.join(datadir, 'shake_data.hdf')
            bu_ic = ShakeMapInputContainer.load(bu_file)
            history = bu_ic.getVersionHistory()
            bu_ic.close()
            if 'history' in history:
                new_line = [timestamp, originator, history['history'][-1][2]]
                history['history'][-1] = new_line
            else:
                history = {'history': []}
                new_line = [timestamp, originator, 1]
                history['history'].append(new_line)
        else:
            #
            # No backup and no existing file. Make this version 1
            #
            history = {'history': []}
            new_line = [timestamp, originator, 1]
            history['history'].append(new_line)

        hdf_file = os.path.join(datadir, 'shake_data.hdf')

        self.logger.debug('Creating input container...')
        shake_data = ShakeMapInputContainer.createFromInput(
            hdf_file,
            config,
            eventxml,
            history,
            rupturefile=rupturefile,
            sourcefile=sourcefile,
            momentfile=momentfile,
            datafiles=datafiles)
        self.logger.debug('Created HDF5 input container in %s' %
                          shake_data.getFileName())
        shake_data.close()
