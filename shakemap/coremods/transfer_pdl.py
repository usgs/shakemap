# stdlib imports
import argparse
import inspect
import os.path
import shutil
import logging

# third party imports
from impactutils.transfer.pdlsender import PDLSender

# local imports
from .transfer_base import TransferBaseModule


class PDLTransfer(TransferBaseModule):
    """
    transfer_pdl - Transfer content via PDL to a remote server.
    """
    command_name = 'transfer_pdl'
    dependencies = [('products/*', False)]

    def __init__(self, eventid):
        # call the parent constructor
        super(PDLTransfer, self).__init__(eventid)

    def execute(self):
        # call parent execute() method
        # this will set the self.info and self.config
        # dictionaries, and self.datadir
        super(PDLTransfer, self).execute()

        # check to see if PDL is a configured method
        if 'pdl' not in self.config:
            logging.info('No PDL transfer has been configured. Returning.')
            return

        # do PDL specific stuff

        pdl_dir = os.path.join(self.datadir, 'pdl')
        products_dir = os.path.join(self.datadir, 'products')
        if not os.path.isdir(pdl_dir):
            raise NotADirectoryError('%s does not exist.' % pdl_dir)

        # get the properties needed for the sender
        properties, product_properties = self.getProperties(self.info)

        downloads_dir = os.path.join(pdl_dir, 'download')
        if os.path.isdir(downloads_dir):
            shutil.rmtree(downloads_dir, ignore_errors=True)
        shutil.copytree(products_dir, downloads_dir)

        # loop over all possible pdl destinations, send products to
        # each one
        for destination, params in self.config['pdl'].items():
            params.update(properties)
            if self.usedevconfig is True:
                if (params['devconfig'] is None or not
                        os.path.isfile(params['devconfig'])):
                    raise FileNotFoundError('Dev config file "%s" does not '
                                            'exist' % params['devconfig'])
                # Swap the config file for the devconfig file
                params['configfile'] = params['devconfig']
                fmt = 'Doing PDL transfer to %s DEV...' % destination
                logging.debug(fmt)
            else:
                fmt = 'Doing PDL transfer to %s...' % destination
                logging.debug(fmt)

            sender = PDLSender(properties=params,
                               local_directory=pdl_dir,
                               product_properties=product_properties)
            if self.cancel:
                msg = sender.cancel()
            else:
                nfiles, msg = sender.send()
                fmt = '%i files sent.  Message from sender: \n"%s"'
                tpl = (nfiles, msg)
                logging.info(fmt % tpl)

        if not self.cancel:
            shutil.rmtree(downloads_dir, ignore_errors=True)

    def parseArgs(self, arglist):
        """
        Set up the object to accept the --dev flag.
        """
        parser = argparse.ArgumentParser(
            prog=self.__class__.command_name,
            description=inspect.getdoc(self.__class__))
        helpstr = ('Send products to the PDL server configured in "devconfig" '
                   'in the transfer.conf configuration file rather than the '
                   'default "configfile".')
        parser.add_argument('-d', '--dev', help=helpstr,
                            action='store_true', default=False)
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
        self.usedevconfig = args.dev
        return args.rem
