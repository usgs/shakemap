import os.path
import tempfile
import textwrap

from configobj import ConfigObj

from shakemap.utils.conf import validate
from shakemap.utils.conf import whatIs
from shakemap.utils.conf import create_default_config

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..'))


def test_conf():
    dummydir = os.path.expanduser('~')
    dummyfile = os.path.join(os.path.expanduser('~'), '.bash_profile')
    data = '''[grind]
    smVs30default = 686.0
    use_gmpe_sc = False
    basin_module = Field2000
    x_grid_interval = 0.025
    y_grid_interval = 0.025
    lonspan = 2.5
    gmroi = 10.0
    iroi = 10.0
    gmdecay = 0.5
    [[gmpe]]
         Zhao06_surface = 0.0,9.9,0,60
         Zhao06_intraslab = 0.0,9.9,60,999
    [[ipe]]
         AW07_CEUS = 0.0,9.9,0,999
    idecay = 0.5
    outlier_deviation_level = 3
    outlier_max_mag = 7.0
    bias_norm = l1
    bias_max_range = 120.0
    bias_min_stations = 6
    bias_max_mag = 7.0
    bias_max_bias = 2.0
    bias_min_bias = -2.0
    bias_log_amp = True
    direct_patch_size = 1000.0
    mi2pgm = WGRW11
    pgm2mi = WGRW11
    [transfer_ftp]
      [[dest1]]
        category = 'webcopy'
        files = *.jpg,
        sendDone = false
        username = user
        password = thispass
        remotehost = ftp.ftptest.org
        remotedirectory = /pub/shakemap
    [transfer_copy]
      [[dest2]]
        category = 'webcopy'
        files = *.jpg,
        sendDone = false
        directory = %s
    [transfer_rsync]
      [[dest2]]
        category = 'webcopy'
        files = *.jpg,
        sendDone = true
        username = fred
        password = password
        privatekey = %s
        remotehost = remotehost.org
        remotedirectory = /home/user/data
    [transfer_pdl]
      [[dest2]]
        category = 'webcopy'
        files = *.jpg,
        sendDone = true
        java = /usr/bin/java
        configfile = %s
        client = %s
        productsource = <SHAKEMAP_NETWORK>
        producttype = shakemap
        productcode = <EVENT_CODE>
        eventsource = <EVENT_NETWORK>
        eventsourcecode = <EVENT_ID>
        privatekey = %s
    ''' % (dummydir, dummyfile, dummyfile, dummyfile, dummyfile)
    fh, tfile = tempfile.mkstemp()
    os.close(fh)
    f = open(tfile, 'wt')
    f.write(textwrap.dedent(data))
    f.close()
    configspec = os.path.join(shakedir, '../shakemap/utils/configspec.ini')
    macros = {'shakemap_network': 'us',
              'event_code': '2015abcd',
              'event_network': 'us',
              'event_id': 'us2015abcd'}
    ret = validate(configspec, tfile, macros=macros)
    config = ConfigObj(tfile)
    configfile = create_default_config(configspec)
    validate(configspec, configfile)
    whatIs(configspec, 'pgm2mi')
    os.remove(tfile)

if __name__ == '__main__':
    test_conf()
    sys.exit(0)
    # where is this script?
    homedir = os.path.dirname(os.path.abspath(__file__))
    configspec = os.path.join(homedir, 'configspec.ini')
    configfile = create_default_config(configspec)

    print(whatIs(configspec, 'pgm2mi'))

    print(validate(configspec, configfile))
