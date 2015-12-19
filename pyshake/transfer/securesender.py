#!/usr/bin/env python

#pip install git+git://github.com/jbardin/scp.py.git
from paramiko import SSHClient
from scp import SCPClient
from sender import Sender,SenderError

class SecureSender(Sender):
    required_props1 = ['privatekey','remotehost','remotedirectory']
    required_props2 = ['username','password','remotehost','remotedirectory']

    def connect(self):
        usePrivateKey = True
        for prop in required_props1:
            if prop not in self.properties.keys():
                usePrivateKey = False
                break
        usePassword = True
        for prop in required_props2:
            if prop not in self.properties.keys():
                usePassword = False
                break
        if not usePrivateKey and not usePassword:
            raise SenderError('Either username/password must be specified, or the name of an SSH private key file.')
        
        ssh = SSHClient()
        if usePrivateKey:
            try:
                ssh.connect(self.properties['remotehost'],
                            key_filename=self.properties['privatekey'],compress=True)
            except Exception,obj:
                raise SenderError('Could not connect with private key file %s' % self.properties['privatekey'])
        else:
            try:
                ssh.connect(self.properties['remotehost'],
                            username=self.properties['username'],password=self.properties['password'],
                            compress=True)
            except Exception,obj:
                raise SenderError('Could not connect with private key file %s' % self.properties['privatekey'])
        return ssh

    def send(self):
        ssh = self.connect()        
        # SCPCLient takes a paramiko transport as its only argument
        scp = SCPClient(ssh.get_transport())
        
        scp.put(self.files, remote_path=self.properties['remotedirectory'])
        scp.put(self.directories, remote_path=self.properties['remotedirectory'],recursive=True)
        nfiles += len(self.files)
        for folder in self.directories:
            nfiles += len(os.listdir(folder))
        scp.close()
        ssh.close()
        return nfiles

    def delete(self):
        ssh = self.connect()
        ssh.exec_command('cd %s;rm %s' % (self.properties['remotedirectory'],' '.join(self.files)))
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
    
if __name__ == '__main__':
    remotehost = sys.argv[1]
    remotefolder = sys.argv[2]
    username = sys.argv[3]
    password = sys.argv[4]
    _testPassword(remotehost,remotefolder,username,password)
    
    
