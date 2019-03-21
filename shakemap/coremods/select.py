"""
Parse STREC output and create/select a GMPE set for an event.
"""
# stdlib imports
import os.path
import shutil
from collections import OrderedDict
import textwrap
import json

# third party imports
from configobj import ConfigObj

# local imports
from .base import CoreModule
import shakemap.utils.config as cfg
from shakemap.utils.probs import get_weights
from shakemap.utils.layers import (validate_config,
                                   update_config_regions)
from shakelib.rupture.origin import Origin


class SelectModule(CoreModule):
    """
    select - Parse STREC output, make a GMPE set, create model_select.conf.
    """

    command_name = 'select'
    targets = [r'model_select\.conf']
    dependencies = [('event.xml', True), ('source.text', False)]
    configs = ['select.conf']

    def execute(self):
        '''
        Parses the output of STREC in accordance with the
        configuration file, creates a new GMPE set specific to the event,
        and writes model_select.conf in the event's 'current' directory.

        Configuration file: select.conf

        Raises:
            NotADirectoryError -- the event's current directory doesn't exist
            FileNotFoundError -- the event.xml file doesn't exist
            ValidateError -- problems with the configuration file
            RuntimeError -- various problems matching the event to a gmpe set
        '''

        # ---------------------------------------------------------------------
        # Get the install and data paths and verify that the even directory
        # exists
        # ---------------------------------------------------------------------
        install_path, data_path = cfg.get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory' % datadir)
        # ---------------------------------------------------------------------
        # Open event.xml and make an Origin object
        # ---------------------------------------------------------------------
        eventxml = os.path.join(datadir, 'event.xml')
        if not os.path.isfile(eventxml):
            raise FileNotFoundError('%s does not exist.' % eventxml)
        momentfile = os.path.join(datadir, 'moment.xml')
        if not os.path.isfile(momentfile):
            momentfile = None
        sourcefile = os.path.join(datadir, 'source.txt')
        if not os.path.isfile(sourcefile):
            sourcefile = None

        org = Origin.fromFile(
            eventxml, sourcefile=sourcefile, momentfile=momentfile)

        #
        # Clear away results from previous runs
        #
        products_path = os.path.join(datadir, 'products')
        if os.path.isdir(products_path):
            shutil.rmtree(products_path, ignore_errors=True)

        # ---------------------------------------------------------------------
        # Get config file from install_path/config, parse and
        # validate it
        # ---------------------------------------------------------------------
        config = ConfigObj(os.path.join(install_path, 'config', 'select.conf'))
        global_data_path = os.path.join(os.path.expanduser('~'),
                                        'shakemap_data')
        validate_config(config, install_path, data_path, global_data_path)

        # ---------------------------------------------------------------------
        # Search through all custom regions, and the first one that we are
        # inside of, take its tectonic region config stuff and replace the
        # default tectonic regions.
        # ---------------------------------------------------------------------
        config = update_config_regions(org.lat, org.lon, config)

        # ---------------------------------------------------------------------
        # Get the default weighting for this event
        # ---------------------------------------------------------------------
        gmmdict, strec_results = get_weights(org, config)

        # ---------------------------------------------------------------------
        # Write the strec results out to a JSON file.
        # ---------------------------------------------------------------------
        jsonfile = os.path.join(datadir, 'strec_results.json')
        with open(jsonfile, 'wt') as f:
            json.dump(strec_results.to_dict(), f)

        # ---------------------------------------------------------------------
        # Create ConfigObj object for output to model_select.conf
        # ---------------------------------------------------------------------
        zc_file = os.path.join(datadir, 'model_select.conf')
        zc_conf = ConfigObj(indent_type='    ')
        zc_conf.filename = zc_file
        zc_conf.initial_comment = textwrap.wrap(
            "This file (model_select.conf) is generated automatically by the "
            "'select' coremod. It will be completely overwritten the next "
            "time select is run. To preserve these settings, or to modify "
            "them, copy this file to a file called 'model.conf' in the "
            "event's current directory. That event-specific model.conf will "
            "be used and model_select.conf will be ignored. (To avoid "
            "confusion, you should probably delete this comment section "
            "from your event-specific model.conf.)", width=75)
        #
        # Add the new gmpe set to the object
        #
        gmpe_set = 'gmpe_' + str(self._eventid) + '_custom'
        zc_conf['gmpe_sets'] = OrderedDict([
            (gmpe_set, OrderedDict([
                ('gmpes', list(gmmdict['gmpelist'])),
                ('weights', list(gmmdict['weightlist'])),
                ('weights_large_dist', 'None'),
                ('dist_cutoff', 'nan'),
                ('site_gmpes', 'None'),
                ('weights_site_gmpes', 'None')
            ]))
        ])
        #
        # Set gmpe to use the new gmpe set
        #
        zc_conf['modeling'] = OrderedDict([
            ('gmpe', gmpe_set),
            ('mechanism', strec_results['FocalMechanism'])
        ])
        if gmmdict['ipe']:
            zc_conf['modeling']['ipe'] = gmmdict['ipe']
        if gmmdict['gmice']:
            zc_conf['modeling']['gmice'] = gmmdict['gmice']
        if gmmdict['ccf']:
            zc_conf['modeling']['ccf'] = gmmdict['ccf']

        zc_conf.write()
