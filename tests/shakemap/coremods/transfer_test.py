#!/usr/bin/env python

import sys
import tempfile
import os.path
import shutil
import argparse

# from shakemap.coremods.transfer import _transfer
from impactutils.io.smcontainers import ShakeMapOutputContainer


def dummy_pdl_test(java, jarfile, keyfile, configfile):
    homedir = os.path.dirname(os.path.abspath(__file__))
    cfile = os.path.join(homedir, '..', '..', 'data', 'containers',
                         'northridge', 'shake_result.hdf')
    products_dir = os.path.join(homedir, '..', '..', 'data',
                                'eventdata', 'northridge', 'current')
    pdl_dir = os.path.join(homedir, '..', '..', 'data', 'eventdata',
                           'northridge', 'current', 'pdl')

    container = ShakeMapOutputContainer.load(cfile)
    config = {'pdl': {'dest1': {'java': java,
                                'jarfile': jarfile,
                                'privatekey': keyfile,
                                'configfile': configfile,
                                'source': 'us'}}}
    # transfermod = TransferModule(eventid)
    _transfer(config, container.getMetadata(), pdl_dir, products_dir)


def dummy_scp_test(remote_host, remote_directory, private_key):
    homedir = os.path.dirname(os.path.abspath(__file__))
    cfile = os.path.join(homedir, '..', '..', 'data', 'containers',
                         'northridge', 'shake_result.hdf')
    products_dir = os.path.join(homedir, '..', '..', 'data',
                                'eventdata', 'northridge', 'current')
    pdl_dir = os.path.join(homedir, '..', '..', 'data', 'eventdata',
                           'northridge', 'current', 'pdl')

    container = ShakeMapOutputContainer.load(cfile)
    config = {'ssh': {'dest1': {'remote_host': remote_host,
                                'remote_directory': remote_directory,
                                'private_key': private_key}}}
    # transfermod = TransferModule(eventid)
    _transfer(config, container.getMetadata(), pdl_dir, products_dir)


def dummy_test_transfer():
    homedir = os.path.dirname(os.path.abspath(__file__))
    cfile = os.path.join(homedir, '..', '..', 'data', 'containers',
                         'northridge', 'shake_result.hdf')
    products_dir = os.path.join(homedir, '..', '..', 'data',
                                'eventdata', 'northridge', 'current')
    pdl_dir = os.path.join(homedir, '..', '..', 'data', 'eventdata',
                           'northridge', 'current', 'pdl')
    eventid = 'ci3144585'

    container = ShakeMapOutputContainer.load(cfile)
    try:
        tdir = tempfile.mkdtemp()
        remote_dir = os.path.join(tdir, eventid)
        config = {'copy': {'local': {'remote_directory': tdir}}}
        # transfermod = TransferModule(eventid)
        _transfer(config, container.getMetadata(), pdl_dir, products_dir)
        nfiles = len(os.listdir(remote_dir))
        nsrcfiles = len(os.listdir(products_dir))
        assert nfiles == nsrcfiles
    except Exception as e:
        print('Exception: %s' % str(e))
    finally:
        if os.path.isdir(tdir):
            shutil.rmtree(tdir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test transfer methods.')
    parser.add_argument(
        '--method',  choices=('pdl', 'scp', 'None'),
        default='None',
        help='Transfer method to test. Default does unit tests')

    # pdl options
    parser.add_argument(
        '--java',
        help='Path to java binary on this system (with pdl method).')
    parser.add_argument(
        '--jarfile',
        help='Path to PDL jar file on this system (with pdl method).')
    parser.add_argument(
        '--keyfile',
        help='Path to PDL private key on this system (with pdl method).')
    parser.add_argument(
        '--configfile',
        help='Path to PDL config file on this system (with pdl method).')

    # scp options
    parser.add_argument(
        '--remote-host',
        help='Remote host (with scp method).')
    parser.add_argument(
        '--remote-directory',
        help='Remote directory (with scp method).')
    parser.add_argument(
        '--private-key',
        help='Path to local private SSH key file (with scp method).')
    args = parser.parse_args()

    if args.method == 'None':
        test_transfer()
    elif args.method == 'pdl':
        if args.java is None or args.jarfile is None or \
                args.keyfile is None or args.configfile is None:
            print('--java, --jarfile, --keyfile, --configfile options must be '
                  'supplied for pdl.')
            sys.exit(1)
        pdl_test(args.java, args.jarfile, args.keyfile, args.configfile)
    elif args.method == 'scp':
        if args.remote_host is None or args.remote_directory is None or \
                args.private_key is None:
            print('--remote-host, --remote-directory, --private-key options '
                  'must be supplied for scp.')
            sys.exit(1)
        scp_test(args.remote_host, args.remote_directory, args.private_key)
