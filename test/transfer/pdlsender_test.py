#!/usr/bin/env python

#stdlib imports
import os.path
import sys
import tempfile
import urllib.request, urllib.error, urllib.parse
import zipfile
import shutil
import io
from distutils import spawn
import textwrap

#third party
from Crypto.PublicKey import RSA

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

#local imports
from shakemap.transfer.pdlsender import PDLSender

def _test_send(internalhub):
    CONFIG = '''senders = sender1
    logdirectory = [FOLDER]/log
    loglevel = FINE
    redirectconsole = false
    enableTracker = false

    [sender1]
    type = gov.usgs.earthquake.distribution.SocketProductSender
    host = [PDLHUB]
    port = 11235'''
    
    PDLURL = 'http://ehppdl1.cr.usgs.gov/ProductClient.zip'
    javabin = spawn.find_executable('java')
    tempdir = None
    try:
        tempdir = tempfile.mkdtemp()
        fh = urllib.request.urlopen(PDLURL)
        zipdata = fh.read()
        fh.close()
        zipf = io.BytesIO(zipdata)
        myzip = zipfile.ZipFile(zipf,'r')
        jarfile = myzip.extract('ProductClient/ProductClient.jar',tempdir)
        myzip.close()
        zipf.close()
        configtext = CONFIG.replace('[PDLHUB]',internalhub)
        configtext = configtext.replace('[FOLDER]',tempdir)
        configfile = os.path.join(tempdir,'config.ini')
        f = open(configfile,'wt')
        f.write(textwrap.dedent(configtext))
        f.close()
        key = RSA.generate(2048)
        keyfile = os.path.join(tempdir,'pdlkey') 
        f = open(keyfile,'wb')
        f.write(key.exportKey('PEM'))
        f.close()
        props = {'java':javabin,
             'jarfile':jarfile,
             'keyfile':keyfile,
             'configfile':configfile,
             'productsource':'ci',
             'producttype':'dummy',
             'productcode':'ci2015abcd',
             'eventsource':'us',
             'eventsourcecode':'us1234abcd'}
        thisfile = os.path.abspath(__file__)
        pdl = PDLSender(properties=props,files=[thisfile])
        pdl.send()
        pdl.delete()
    except Exception as obj:
        pass
    #remove temporary pdl folder with jarfile, config, and keyfile in it
    if tempdir is not None:
        shutil.rmtree(tempdir)
    
if __name__ == '__main__':
    internalhub = sys.argv[1] #this needs to be the hostname of a PDL server that does not require a registered public key
    _test_send(internalhub)
