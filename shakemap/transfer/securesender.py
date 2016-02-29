#!/usr/bin/env python

#stdlib imports
import sys
import os.path

#depends on https://github.com/jbardin/scp.py
#pip install git+git://github.com/jbardin/scp.py.git
from paramiko import SSHClient
from scp import SCPClient
from sender import Sender

#local imports
from shakemap.utils.exception import ShakeMapException

class SecureSender(Sender):
    '''Class for sending and deleting files and directories via SSH.
    '''
    required_props1 = ['privatekey','remotehost','remotedirectory']
    required_props2 = ['username','password','remotehost','remotedirectory']

    def connect(self):
        """Initiate an ssh connection with properties passed to constructor.
        :returns:
          Instance of the paramiko SSHClient class.
        """
        usePrivateKey = True
        for prop in self.required_props1:
            if prop not in list(self.properties.keys()):
                usePrivateKey = False
                break
        usePassword = True
        for prop in self.required_props2:
            if prop not in list(self.properties.keys()):
                usePassword = False
                break
        if not usePrivateKey and not usePassword:
            raise ShakeMapException('Either username/password must be specified, or the name of an SSH private key file.')
        
        ssh = SSHClient()
        #load hosts found in ~/.ssh/known_hosts
        ssh.load_system_host_keys() #should we not assume that the user has these configured already?
        if usePrivateKey:
            try:
                ssh.connect(self.properties['remotehost'],
                            key_filename=self.properties['privatekey'],compress=True)
            except Exception as obj:
                raise ShakeMapException('Could not connect with private key file %s' % self.properties['privatekey'])
        else:
            try:
                ssh.connect(self.properties['remotehost'],
                            username=self.properties['username'],password=self.properties['password'],
                            compress=True)
            except Exception as obj:
                raise ShakeMapException('Could not connect with private key file %s' % self.properties['privatekey'])
        return ssh

    def send(self):
        '''Send any files or folders that have been passed to constructor.
        :returns:
          Number of files sent to remote FTP server.
        '''
        nfiles = 0
        ssh = self.connect()        
        # SCPCLient takes a paramiko transport as its only argument
        scp = SCPClient(ssh.get_transport())
        try:
            if len(self.files):
                scp.put(self.files, remote_path=self.properties['remotedirectory'])
            if len(self.directories):
                scp.put(self.directories, remote_path=self.properties['remotedirectory'],recursive=True)
        except Exception as obj:
            pass
        nfiles += len(self.files)
        for folder in self.directories:
            nfiles += len(os.listdir(folder))
        scp.close()
        ssh.close()
        return nfiles

    def delete(self):
        '''Delete any files and folders that have been passed to constructor.
        :returns:
          The number of files deleted on remote SSH server.
        '''
        ssh = self.connect()
        rfiles = []
        for filename in self.files:
            fbase,fname = os.path.split(filename)
            rdic = self.properties['remotedirectory']
            rfile = os.path.join(rdic,fname)
            rfiles.append(rfile)
        cmd = 'rm %s' % (' '.join(rfiles))
        ssh.exec_command(cmd)
        ssh.close()


def _testPassword(remotehost,remotefolder,username,password):
    props = {'remotehost':remotehost,
             'remotedirectory':remotefolder,
             'username':username,
             'password':password}
    thisfile = os.path.abspath(__file__)
    securesend = SecureSender(properties=props,files=[thisfile])
    securesend.send()
    securesend.delete()

def _testKey(remotehost,remotefolder,privatekey):
    props = {'remotehost':remotehost,
             'remotedirectory':remotefolder,
             'privatekey':privatekey}
    thisfile = os.path.abspath(__file__)
    securesend = SecureSender(properties=props,files=[thisfile])
    securesend.send()
    securesend.delete()
    
if __name__ == '__main__':
    if len(sys.argv) == 5:
        remotehost = sys.argv[1]
        remotefolder = sys.argv[2]
        username = sys.argv[3]
        password = sys.argv[4]
        _testPassword(remotehost,remotefolder,username,password)
    else:
        remotehost = sys.argv[1]
        remotefolder = sys.argv[2]
        privatekey = sys.argv[3]
        _testKey(remotehost,remotefolder,privatekey)
    
    
    
    
