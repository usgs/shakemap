import os.path
from configobj import ConfigObj

def get_data_path():
    homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
    datadir = os.path.join(homedir,'..','data')
    return datadir

def get_config_paths():
    config_file = os.path.join(os.path.expanduser('~'),'.shakemap','profiles.conf')
    config = ConfigObj(config_file)
    profile_name = config['profile']
    profile = config['profiles'][profile_name]
    install = profile['install_path']
    data = profile['data_path']
    return (install,data)
