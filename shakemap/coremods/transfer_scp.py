# stdlib imports
import os.path
from datetime import datetime
import shutil
import logging

# third party imports
from impactutils.transfer.securesender import SecureSender

# local imports
from .transfer_base import TransferBaseModule


class SCPTransfer(TransferBaseModule):
    """
    transfer_pdl - Transfer content via SCP to a remote server.
    """
    command_name = 'transfer_scp'
    dependencies = [('products/*', False)]

    def __init__(self, eventid):
        # call the parent constructor
        super(SCPTransfer, self).__init__(eventid)

    def execute(self):
        # call parent execute() method
        # this will set the self.info and self.config
        # dictionaries, and self.datadir
        super(SCPTransfer, self).execute()

        # check to see if SCP is a configured method
        if 'scp' not in self.config:
            logging.info('No SCP transfer has been configured. Returning.')
            return

        # get the properties needed for the sender
        properties, product_properties = self.getProperties(self.info)

        # get the products directory
        product_dir = os.path.join(self.datadir, 'products')

        # loop over all possible scp destinations, send products to
        # each one
        for destination, params in self.config['scp'].items():
            # append the event ID to the remote_directory
            pdir = params['remote_directory']
            params['remote_directory'] = os.path.join(pdir, self._eventid)

            params.update(properties)
            fmt = 'Doing SCP transfer to %s...' % destination
            logging.debug(fmt)

            sender = SecureSender(properties=params,
                                  local_directory=product_dir)
            if self.cancel:
                msg = sender.cancel()
            else:
                try:
                    nfiles, msg = sender.send()
                except Exception as e:
                    logging.warning(str(e))
                    raise(e)
                fmt = '%i files sent.  Message from sender: \n"%s"'
                tpl = (nfiles, msg)
                logging.info(fmt % tpl)
