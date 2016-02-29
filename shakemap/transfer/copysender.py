#!/usr/bin/env python

#stdlib imports
import shutil
import os.path
import tempfile

#local imports
from sender import Sender

#local imports
from shakemap.utils.exception import ShakeMapException

class CopySender(Sender):
    '''(Really simple) Class for sending and deleting files and directories via system copy and delete.
    '''
    def send(self):
        '''Send any files or folders that have been passed to constructor.
        :returns:
          Number of files sent to local directory.
        '''
        nfiles = 0
        if 'directory' not in self.properties:
            raise ShakeMapException('Property "directory" not specified.')
        
        if not os.path.isdir(self.properties['directory']):
            raise ShakeMapException('Output directory "%s" does not exist.' % self.properties['directory'])
        
        for filename in self.files:
            shutil.copy(filename,self.properties['directory'])
            
        nfiles += len(self.files)
        
        for folder in self.directories:
            shutil.copytree(folder,self.properties['directory'])
            nfiles += len(os.walk(folder).next()[2])

        return nfiles

    def delete(self):
        '''Delete any files and folders that have been passed to constructor.
        :returns:
          The number of files deleted from local directory.
        '''
        if 'directory' not in self.properties:
            raise ShakeMapException('Property "directory" not specified.')
        
        if not os.path.isdir(self.properties['directory']):
            raise ShakeMapException('Output directory "%s" does not exist.' % self.properties['directory'])

        for filename in self.files:
            fbase,fname = os.path.split(filename)
            dfile = os.path.join(self.properties['directory'],fname)
            os.remove(dfile)

        for folder in self.directories:
            fbase,dirname = os.path.split(folder)
            dfolder = os.path.join(self.properties['directory'],dirname)
            shutil.rmtree(dfolder)


def _test():
    thisfile = os.path.abspath(__file__)
    tempdir = tempfile.mkdtemp()
    try:
        cpsender = CopySender(properties={'directory':tempdir},files=[thisfile])
        nfiles = cpsender.send()
        nfiles = cpsender.delete()
    except Exception as obj:
        raise ShakeMapException('Failed to copy or delete a file.')
    shutil.rmtree(tempdir)

if __name__ == '__main__':
    _test()
        
        
    
