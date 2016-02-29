#!/usr/bin/env python

#stdlib imports
import os.path
import tempfile
import subprocess
import sys
import urllib.request, urllib.error, urllib.parse
import zipfile
import shutil
import io
from distutils import spawn
import textwrap

#third party
from Crypto.PublicKey import RSA
from neicio.cmdoutput import getCommandOutput

#local
from sender import Sender

#local imports
from shakemap.utils.exception import ShakeMapException

class PDLSender(Sender):
    pdlcmd = '[JAVA] -jar [JARFILE] --send --status=[STATUS] --source=[PRODUCTSOURCE] --type=[PRODUCTTYPE] --code=[PRODUCTCODE] --eventsource=[EVENTSOURCE] --eventsourcecode=[EVENTSOURCECODE] --privateKey=[KEYFILE]  --configFile=[CONFIGFILE] [FILE] [DIRECTORY]'
    required_properties = ['java','jarfile','keyfile','configfile',
                           'productsource','producttype','productcode',
                           'eventsource','eventsourcecode']
    def delete(self):
        for prop in self.required_properties:
            if prop not in list(self.properties.keys()):
                raise ShakeMapException('"%s" property must be supplied to send via PDL')

        #build pdl command line from properties
        self.properties['status'] = 'DELETE'
        self.properties['files'] = ''
        self.properties['directories'] = ''
        cmd = self.pdlcmd
        for propkey,propvalue in self.properties.items():
            cmd = cmd.replace('['+propkey.upper()+']',propvalue)
        
        retcode,stdout,stderr = getCommandOutput(cmd)
        if not retcode:
            fmt = 'Could not delete product "%s" due to error "%s"'
            tpl = (code,stdout+stderr)
            raise ShakeMapException(fmt % tpl)
        
    def send(self):
        #we can really only support sending of one file and/or one directory, so error out
        #if someone has specified more than one of either.
        if len(self.files) > 1:
            raise ShakeMapException('For PDL, you may only send one file at a time.')
        if len(self.directories) > 1:
            raise ShakeMapException('For PDL, you may only send one directory at a time.')

        #make sure we have all the required properties
        for prop in self.required_properties:
            if prop not in list(self.properties.keys()):
                raise ShakeMapException('"%s" property must be supplied to send via PDL')

        
        #build pdl command line from properties
        self.properties['command'] = 'send'
        self.properties['status'] = 'UPDATE'
        if self.files:
            self.properties['file'] = self.files[0]
        else:
            self.properties['file'] = ''
        if self.directories:
            self.properties['directory'] = self.directories[0]
        else:
            self.properties['directory'] = ''
        cmd = self.pdlcmd
        for propkey,propvalue in self.properties.items():
            cmd = cmd.replace('['+propkey.upper()+']',propvalue)

        #call PDL on the command line
        retcode,stdout,stderr = getCommandOutput(cmd)
        if not retcode:
            fmt = 'Could not send product "%s" due to error "%s"'
            tpl = (code,stdout+stderr)
            raise ShakeMapException(fmt % tpl)

        #return the number of files we just sent
        nfiles = 0
        if self.properties['file']:
            nfiles += 1
        if self.properties['directory']:
            nfiles += len(os.listdir(self.properties['directory']))

        return nfiles
                        

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
        zipf = io.StringIO(zipdata)
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
        f = open(keyfile,'wt')
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
