# stdlib imports
import os.path
from datetime import datetime
import shutil
import logging

# third party imports
from impactutils.transfer.emailsender import EmailSender
import configobj

# local imports
from .transfer_base import TransferBaseModule
from shakemap.utils.macros import get_macros


class EmailTransfer(TransferBaseModule):
    """
    transfer_email - Transfer content via Email.
    """
    command_name = 'transfer_email'
    dependencies = [('products/*', False)]

    def __init__(self, eventid):
        # call the parent constructor
        super(EmailTransfer, self).__init__(eventid)

    def execute(self):
        # call parent execute() method
        # this will set the self.info and self.config
        # dictionaries, and self.datadir
        super(EmailTransfer, self).execute()

        # check to see if email is a configured method
        if 'email' not in self.config:
            logging.info('No email transfer has been configured. Returning.')
            return

        # see what the user set for the mail_once setting for all destinations
        mail_once = self.config['email']['mail_once']

        # then check for the mail flag file
        mailfile = os.path.join(self.datadir, '.mailed')
        if os.path.isfile(mailfile) and not self.cancel:
            msg = 'Mail has already been generated for this event. Returning.'
            logging.info(msg)
            return

        # get the properties needed for the sender
        properties, product_properties = self.getProperties(self.info)

        # get the products directory path
        product_dir = os.path.join(self.datadir, 'products')

        # get the macros that may be in the email sender config
        macros = get_macros(self.info)

        # loop over all possible email destinations, send products to
        # each one
        for destination, params in self.config['email'].items():
            if not isinstance(params, configobj.Section):
                continue
            params.update(properties)
            fmt = 'Doing email transfer to %s...' % destination
            logging.debug(fmt)

            # replace macro strings with actual strings
            for pkey, param in params.items():
                for macro, replacement in macros.items():
                    if isinstance(param, str):
                        try:
                            param = param.replace('[%s]' % macro,
                                                  replacement)
                            params[pkey] = param
                        except Exception as e:
                            x = 1

            # get full path to all file attachments
            attachments = []
            for lfile in params['attachments']:
                fullfile = os.path.join(product_dir, lfile)
                if not os.path.isfile(fullfile):
                    logging.warn('%s does not exist.' % fullfile)
                attachments.append(fullfile)

            sender = EmailSender(properties=params,
                                 local_files=attachments)
            if self.cancel:
                msg = sender.cancel()
            else:
                try:
                    nfiles, msg = sender.send()
                    with open(mailfile, 'wt') as f:
                        f.write('Mailed at %s UTC' % (str(datetime.utcnow())))
                except Exception as e:
                    logging.warning(str(e))
                    raise(e)
                fmt = '%i files sent.  Message from sender: \n"%s"'
                tpl = (nfiles, msg)
                logging.info(fmt % tpl)
