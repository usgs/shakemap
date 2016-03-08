#!/usr/bin/env python

#stdlib imports
import os.path
import sys
import tempfile
import shutil

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

from shakemap.transfer.copysender import CopySender

def _test():
    print('Testing basic file system copy...')
    thisfile = os.path.abspath(__file__)
    tempdir = tempfile.mkdtemp()
    try:
        cpsender = CopySender(properties={'directory':tempdir},files=[thisfile])
        nfiles = cpsender.send()
        nfiles = cpsender.delete()
    except Exception as obj:
        raise SenderError('Failed to copy or delete a file.')
    shutil.rmtree(tempdir)
    print('Passed basic file system copy.')

if __name__ == '__main__':
    _test()
        
        
    
