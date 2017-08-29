#!/usr/bin/env python

import argparse
import os.path
import sys
import pathlib
import glob
import re

from impactutils.io.cmd import get_command_output

def main(args):

    #-------------------------------------------------------------
    # Some additional useful directories
    #-------------------------------------------------------------
    REPO_DIR = os.path.dirname(os.path.abspath(__file__))
    PACKAGE_DIR = os.path.join(REPO_DIR, 'shakemap')
    DOC_SRC_DIR = os.path.join(REPO_DIR, 'doc_source')
    API_DIR = os.path.join(DOC_SRC_DIR, 'apidoc')
    DOCS_DIR = os.path.join(REPO_DIR, 'docs')

    #-------------------------------------------------------------
    # get the human-friendly version of the ShakeMap version
    #-------------------------------------------------------------
    verstr = '4.0a'

    #-------------------------------------------------------------
    # what is the package called and who are the authors
    #-------------------------------------------------------------
    PACKAGE = "ShakeMap 4.0a API"
    AUTHORS = 'Bruce Worden, Eric Thompson, Mike Hearne'

    #-------------------------------------------------------------
    # run the api doc command; this creates the .rst files
    #-------------------------------------------------------------
    sys.stderr.write('Building shakemap API documentation (REST)...\n')
    sphinx_cmd = 'sphinx-apidoc -o %s -f -e -d 12 -H "%s" -A "%s"'\
                 ' -V %s -T %s' % (API_DIR, PACKAGE, AUTHORS, verstr,
                                PACKAGE_DIR)
    res, stdout, stderr = get_command_output(sphinx_cmd)

    if not res:
        raise Exception('Could not build ShakeMap API documentation'
                        ' - error "%s".' % stderr.decode())
    if args.verbose:
        print(stdout.decode('utf-8'))
        print(stderr.decode('utf-8'))

    #--------------------------------------------
    # try to clean up some of the excess labeling
    #--------------------------------------------
#    clean_cmd = "sed -e 's/ module//g' -i '' `find %s/*.rst -type f "\
#                "-maxdepth 0 -print`" % API_DIR
#    res, stdout, stderr = get_command_output(clean_cmd)
#    clean_cmd = "sed -e 's/ package//g' -i '' `find %s/*.rst -type f "\
#                "-maxdepth 0 -print`" % API_DIR
#    res, stdout, stderr = get_command_output(clean_cmd)
#    clean_cmd = "sed -e '/Subpackages/d' -i '' `find %s/*.rst -type f "\
#                "-maxdepth 0 -print`" % API_DIR
#    res, stdout, stderr = get_command_output(clean_cmd)
#    clean_cmd = "sed -e '/-.*-/d' -i '' `find %s/*.rst -type f "\
#                "-maxdepth 0 -print`" % API_DIR
#    res, stdout, stderr = get_command_output(clean_cmd)

    #-------------------------------------------------------------
    # Go to the api directory and build the html
    #-------------------------------------------------------------
    sys.stderr.write('Building shakemap manual (HTML)...\n')
    res, stdout, stderr = get_command_output('sphinx-build -a -E %s %s'
                                             % (DOC_SRC_DIR, DOCS_DIR))
    if not res:
        raise Exception('Could not build HTML. - '
                        'error "%s"' % stderr.decode())
    if args.verbose:
        print(stdout.decode('utf-8'))
        print(stderr.decode('utf-8'))

    pathlib.Path(os.path.join(DOCS_DIR, '.nojekyll')).touch(exist_ok=True)


if __name__ == '__main__':
    desc = 'Create API documentation for ShakeMap and optionally post '\
           'HTML output to GitHub.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Produce more output to the screen. ')

    pargs = parser.parse_args()
    main(pargs)
