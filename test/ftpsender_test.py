#!/usr/bin/env python

#stdlib imports
import os.path
import sys
import urllib.request, urllib.error, urllib.parse

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

#local imports
from shakemap.transfer.ftpsender import FTPSender
from shakemap.utils.exception import ShakeMapException

def _testSendFile(properties):
    print('Testing sending single file...')
    thisfile = os.path.abspath(__file__)
    thispath,thisfilename = os.path.split(thisfile)
    try:
        sender = FTPSender(properties=properties,files=[thisfile])
        sender.send()
        url = 'ftp://%s%s%s' % (properties['host'],properties['directory'],thisfilename)
        fh = urllib.request.urlopen(url)
        fh.close()
        sender.delete()
    except Exception as obj:
        fmt = 'Test failed - you may have a file called %s on host %s and directory %s'
        tpl = (thisfile,properties['host'],['directory'])
        raise ShakeMapException(fmt % tpl)
    print('Passed sending single file.')

def _testSendFolder(properties):
    #modify this to create a temporary folder and send that - I think __pycache__ is screwing up the deletes...
    #although maybe I should test deleting directories with directories in them...
    print('Testing sending folder...')
    thisfile = os.path.abspath(__file__)
    thispath,thisfilename = os.path.split(thisfile)
    try:
        sender = FTPSender(properties=properties,directory=thispath)
        sender.send()
        url = 'ftp://%s%s' % (properties['host'],properties['directory'])
        fh = urllib.request.urlopen(url)
        fh.close()
        sender.delete()
    except Exception as obj:
        fmt = 'Test failed - you may have a file called %s on host %s and directory %s'
        tpl = (thisfile,properties['host'],['directory'])
        raise ShakeMapException(fmt % tpl)
    print('Passed sending folder.')
    
if __name__ == '__main__':
    #try logging into an FTP server that supports anonymous login
    host = sys.argv[1]
    folder = sys.argv[2]
    user = ''
    password = 'user@anonymous.org'
    props = {'host':host,
             'directory':folder,
             'user':user,
             'password':password}
    _testSendFile(props)
    _testSendFolder(props)
