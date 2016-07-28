#!/usr/bin/env python

# stdlib imports
import sys
import os.path

# hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
# put this at the front of the system path, ignoring any installed mapio stuff
sys.path.insert(0, shakedir)

from shakemap.transfer.securesender import SecureSender


def _testPassword(remotehost, remotefolder, username, password):
    props = {'remotehost': remotehost,
             'remotedirectory': remotefolder,
             'username': username,
             'password': password}
    thisfile = os.path.abspath(__file__)
    securesend = SecureSender(properties=props, files=[thisfile])
    securesend.send()
    securesend.delete()


def _testKey(remotehost, remotefolder, privatekey):
    props = {'remotehost': remotehost,
             'remotedirectory': remotefolder,
             'privatekey': privatekey}
    thisfile = os.path.abspath(__file__)
    securesend = SecureSender(properties=props, files=[thisfile])
    securesend.send()
    securesend.delete()

if __name__ == '__main__':
    if len(sys.argv) == 5:
        remotehost = sys.argv[1]
        remotefolder = sys.argv[2]
        username = sys.argv[3]
        password = sys.argv[4]
        _testPassword(remotehost, remotefolder, username, password)
    else:
        remotehost = sys.argv[1]
        remotefolder = sys.argv[2]
        privatekey = sys.argv[3]
        _testKey(remotehost, remotefolder, privatekey)
