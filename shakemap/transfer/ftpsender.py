#!/usr/bin/env python

#stdlib imports
from ftplib import FTP
import os.path
import sys
import urllib.request, urllib.error, urllib.parse

#local
from sender import Sender,SenderError
    
class FTPSender(Sender):
    '''Class for sending and deleting files and directories via FTP.
    '''
    def setup(self):
        """
        Initiate an ftp connection with properties passed to constructor.
        :returns:
          Instance of the ftplib.FTP class.
        """
        if 'host' not in list(self.properties.keys()):
            raise NameError('"host" keyword must be supplied to send via FTP')
        if 'directory' not in list(self.properties.keys()):
            raise NameError('"directory" keyword must be supplied to send via FTP')
        host = self.properties['host']
        folder = self.properties['directory']
        try:
            dirparts = folder.strip().split('/')
            ftp = FTP(host)
            if 'user' in self.properties:
                user = self.properties['user']
            else:
                user = ''
            if 'password' in self.properties:
                password = self.properties['password']
            else:
                password = ''
            if user == '':
                ftp.login()
            else:
                ftp.login(user,password)
            for d in dirparts:
                if d == '':
                    continue
                try:
                    ftp.cwd(d)
                except ftplib.error_perm as msg:
                    raise SenderError('Could not login to host "%s" and navigate to directory "%s"' % (host,folder))
        except Exception as obj:
            raise SenderError('Could not send to %s.  Error "%s"' % (host,str(obj)))
        return ftp

    def delete(self):
        '''Delete any files and folders that have been passed to constructor.
        :returns:
          The number of files deleted on remote FTP server.
        '''
        ftp = self.setup()
        nfiles = 0
        host = self.properties['host']
        folder = self.properties['directory']
        if self.files is not None:
            for f in self.files:
                fbase,fpath = os.path.split(f)
                ftp.delete(fpath)
                nfiles += 1
        if self.directories is not None:
            for directory in self.directories:
                root,thisfolder = os.path.split(directory) #root is the top level local directory
                for path, subdirs, files in os.walk(directory):
                    mpath = path.replace(root,'').lstrip(os.sep) #mpath is the relative path on the ftp server
                    allfiles = ftp.nlst()
                    if mpath not in allfiles:
                        print('Could not find directory %s on ftp server.' % mpath)
                        continue
                    ftpfolder = os.path.join(folder,mpath) #full path to the folder on ftp server
                    ftp.cwd(ftpfolder)
                    for f in files:
                        #f is the file name within the current folder
                        ftp.delete(f)
                        nfiles += 1
                    ftp.cwd(folder) #go back to the root 
                    ftp.rmd(ftpfolder)
        ftp.quit()
        return nfiles
    
    def send(self):
        '''Send any files or folders that have been passed to constructor.
        :returns:
          Number of files sent to remote SSH server.
        '''
        if 'host' not in list(self.properties.keys()):
            raise NameError('"host" keyword must be supplied to send via FTP')
        if 'directory' not in list(self.properties.keys()):
            raise NameError('"directory" keyword must be supplied to send via FTP')
        try:
            host = self.properties['host']
            folder = self.properties['directory']
            ftp = self.setup()
            #ftp.cwd(self.properties['directory'])
            nfiles = 0
            if self.files is not None:
                for f in self.files:
                    self.__sendfile(f,ftp)
                    nfiles += 1
            if self.directories is not None:
                for directory in self.directories:
                    root,thisfolder = os.path.split(directory) #root is the top level local directory
                    for path, subdirs, files in os.walk(directory):
                        #mpath is the relative path on the ftp server
                        mpath = path.replace(root,'').lstrip(os.sep) 
                        allfiles = ftp.nlst()
                        if mpath not in allfiles:
                            ftp.mkd(mpath)
                        ftpfolder = os.path.join(folder,mpath) #full path to the folder on ftp server
                        ftp.cwd(ftpfolder)
                        for f in files:
                            #f is the file name within the current folder
                            fpath = os.path.join(path,f) #the full path to the local file
                            self.__sendfile(fpath,ftp)
                            nfiles += 1
                        ftp.cwd(folder) #go back to the root 
            ftp.quit()
            return nfiles
                    
        except Exception as obj:
            raise SenderError('Could not send to %s.  Error "%s"' % (host,str(obj)))

    def __sendfile(self,filename,ftp):
        '''Internal function used to send a file using an FTP object.
        :param filename:
          Local filename
        :param ftp:
          Instance of FTP object.
        '''
        fbase,fpath = os.path.split(filename)
        cmd = "STOR " + fpath #we don't tell the ftp server about the local path to the file
        ftp.storbinary(cmd,open(filename,"rb"),1024) #actually send the file

def _testSendFile(properties):
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
        raise SenderError(fmt % tpl)

def _testSendFolder(properties):
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
        raise SenderError(fmt % tpl)

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
