#!/usr/bin/env python

import argparse
import os.path
import shutil
import sys
import re
from distutils.dir_util import copy_tree

from shakemap.utils.misc import get_command_output


DEFAULT_TAG = '0.1'

def main(args):

    #-------------------------------------------------------------
    # where should .rst files, Makefile, _build folder be written?
    #-------------------------------------------------------------
    API_DIR = os.path.join(os.path.expanduser('~'), '__api-doc')
    shutil.rmtree(API_DIR, ignore_errors=True)

    #-------------------------------------------------------------
    # where should the temporary clone of the shakemap gh-pages repo live?
    #-------------------------------------------------------------
    TOP_DIR = os.path.join(os.path.expanduser('~'), '__shake-doc')
    CLONE_DIR = os.path.join(TOP_DIR, 'html')
    shutil.rmtree(TOP_DIR, ignore_errors=True)

    #-------------------------------------------------------------
    # Some additional useful directories
    #-------------------------------------------------------------
    REPO_DIR = os.path.dirname(os.path.abspath(__file__))
    DOC_DIR = os.path.join(REPO_DIR, 'doc')
    REST_DIR = os.path.join(REPO_DIR, 'rest')
    PACKAGE_DIR = os.path.join(REPO_DIR, 'shakemap')

    #-------------------------------------------------------------
    # get the human-friendly version of the ShakeMap version
    #-------------------------------------------------------------
    res, verstr, stderr = get_command_output('git describe --tags')
    verstr = verstr.decode().strip()
    if not len(verstr):
        verstr = DEFAULT_TAG
    else:
        if verstr[0] == 'v':
            verstr = verstr[1:]
        spl = re.findall(r"[\w']+", verstr)
        verstr = "%s.%s.%s" % (spl[0], spl[1], spl[2])

    #-------------------------------------------------------------
    # what is the package called and who are the authors
    #-------------------------------------------------------------
    PACKAGE = "shakemap"
    AUTHORS = 'Bruce Worden, Eric Thompson, Mike Hearne'

    # find the make command on this system
    res, stdout, stderr = get_command_output('which make')
    if not res:
        print('Could not find the "make" command on your system. Exiting.')
        sys.exit(1)
    make_cmd = stdout.decode().strip()

    #-------------------------------------------------------------
    # clone the repository
    #-------------------------------------------------------------
    if args.post:
        sys.stderr.write('Cloning shakemap gh-pages branch...\n')
        if os.path.isdir(CLONE_DIR):
            shutil.rmtree(CLONE_DIR)
        clonecmd = 'git clone -b gh-pages https://github.com/usgs/'\
                   'shakemap.git %s' % CLONE_DIR
        res, stdout, stderr = get_command_output(clonecmd)
        if not res:
            raise Exception('Could not clone gh-pages branch.')

        # Delete everything in the repository (except hidden git files)
        cmd = 'rm -fr %s' %CLONE_DIR

    #-------------------------------------------------------------
    # run the api doc command; this creates the .rst files
    #-------------------------------------------------------------
    sys.stderr.write('Building shakemap API documentation (REST)...\n')
    sphinx_cmd = 'sphinx-apidoc -o %s -f -e -l -d 12 -F -H %s -A "%s"'\
                 ' -V %s %s' % (API_DIR, PACKAGE, AUTHORS, verstr,
                                PACKAGE_DIR)

    res, stdout, stderr = get_command_output(sphinx_cmd)

    if not res:
        raise Exception('Could not build ShakeMap API documentation'
                        ' - error "%s".' % stderr)

    #-------------------------------------------------------------
    # change index.rst to api_index.rst
    #-------------------------------------------------------------
    move_cmd = 'mv %s/index.rst %s/api_index.rst' % (API_DIR, API_DIR)
    res, stdout, stderr = get_command_output(move_cmd)

    # Change name of API documentation in api_index.rst
    cmd = "sed -i -e 's/Welcome to shakemap.*/ShakeMap 4.0 API/g' "\
          "%s/api_index.rst" % API_DIR
    res, stdout, stderr = get_command_output(cmd)

    #--------------------------------------------
    # try to clean up some of the excess labeling
    #--------------------------------------------
    clean_cmd = "sed -i -e 's/ module//g' `find %s/*.rst -type f "\
                "-maxdepth 0 -print`" % API_DIR
    res, stdout, stderr = get_command_output(clean_cmd)
    clean_cmd = "sed -i -e 's/ package//g' `find %s/*.rst -type f "\
                "-maxdepth 0 -print`" % API_DIR
    res, stdout, stderr = get_command_output(clean_cmd)
    clean_cmd = "sed -i -e '/Subpackages/d' `find %s/*.rst -type f "\
                "-maxdepth 0 -print`" % API_DIR
    res, stdout, stderr = get_command_output(clean_cmd)
    clean_cmd = "sed -i -e '/-.*-/d' `find %s/*.rst -type f "\
                "-maxdepth 0 -print`" % API_DIR
    res, stdout, stderr = get_command_output(clean_cmd)

    #-------------------------------------------------------------
    # run the make command to build the shakemap manual (pdf version)
    #-------------------------------------------------------------
    if not args.nopdf:
        sys.stderr.write('Building shakemap manual (PDF)...\n')
        os.chdir(DOC_DIR)
        # Need to run HTML to create __shakedoc/html/_static
        manualcmd = '%s latexpdf' % make_cmd
        res, stdout, stderr = get_command_output(manualcmd)
        if not res:
            raise Exception('Could not build the PDF version of the ShakeMap '
                            'manual - error "%s".' % stderr)

        # Move PDF to API_DIR/_static
        spdf = '%s/latex/ShakeMapManual.pdf' %TOP_DIR
        dpdf = '%s/_static/ShakeMapManual.pdf' %API_DIR
        shutil.copy2(spdf, dpdf)

    #-------------------------------------------------------------
    # Edit the conf.py file to include the theme.
    #-------------------------------------------------------------
    fname = os.path.join(API_DIR, 'conf.py')
    f = open(fname, 'at')
    f.write("sys.path.insert(0, os.path.abspath('%s'))\n" % (REPO_DIR))

    #-------------------------------------
    # Built in theme:
    #-------------------------------------
#    f.write("html_theme = 'haiku'\n")
    #-------------------------------------

    #-------------------------------------
    # RTD theme
    #-------------------------------------
    f.write("import sphinx_rtd_theme\n")
    f.write("html_theme = 'sphinx_rtd_theme'\n")
    f.write("html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]\n")
    f.write("html_theme_options = {\n")
    f.write("    'collapse_navigation': False,\n")
    f.write("}\n")
    #-------------------------------------

    #-------------------------------------
    # Bootstrap theme
    #-------------------------------------
#    f.write("import sphinx_bootstrap_theme\n")
#    f.write("html_theme = 'bootstrap'\n")
#    f.write("html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()\n")
#    f.write("html_theme_options = {\n")
#    f.write("    'bootswatch_theme': \"Sandstone\",\n")
#    f.write("    'navbar_site_name': \"API Documentation\",\n")
#    f.write("    'navbar_sidebarrel': False,\n")
#    f.write("    'navbar_pagenav': False,\n")
#    f.write("    'navbar_links': [\n")
#    f.write("        (\"Manual\", \"../shake_index.html\", 1),\n")
#    f.write("    ],\n")
#    f.write("}\n")
    #-------------------------------------

    # Napolean extension? Supports Goggle and Numpy style docstrings, but it
    # also has some side effects such as restrictions on what sections are
    # allowed and it seems to suppress the [source] link to code; maybe this
    # is a configurable option though.
#    f.write("extensions = ['sphinx.ext.autodoc', 'sphinxcontrib.napoleon']\n")

    # This line is needed to inclue __init__ methods in documentation
    f.write("autoclass_content = 'both'\n")
    f.write("html_show_copyright = False\n")
    f.close()

    #-------------------------------------------------------------
    # Copy the manual REST files to the API directory
    #-------------------------------------------------------------

    # put aside Makefile so it doesn't get overwritten
    oldmake = '%s/Makefile' %API_DIR
    tmpmake = '%s/Makefile_save' %API_DIR
    os.rename(oldmake, tmpmake)

    # move files inot API directory; this should raise exceptions if any
    # files will get overwritten. 
    copy_tree(DOC_DIR, API_DIR)

    # put Makefile back
    os.rename(tmpmake, oldmake)

    # Move index.rst to manual_index.rst
    ind1 = '%s/index.rst' %API_DIR
    ind2 = '%s/manual_index.rst' %API_DIR
    os.rename(ind1, ind2)

    # Copy index.rst from rest/ directory into build directory
    restind = '%s/index.rst' %REST_DIR
    apiind = '%s/index.rst' %API_DIR
    shutil.copy2(restind, apiind)

    #-------------------------------------------------------------
    # Go to the api directory and build the html
    #-------------------------------------------------------------
    sys.stderr.write('Building shakemap manual (HTML)...\n')
    os.chdir(API_DIR)
    res, stdout, stderr = get_command_output('%s html' % make_cmd)
    if not res:
        raise Exception('Could not build HTML for API documentation. - '
                        'error "%s"' % stderr)

    #-------------------------------------------------------------
    # Copy the generated content to the gh-pages branch we created
    # earlier
    #-------------------------------------------------------------
    htmldir = os.path.join(API_DIR, '_build', 'html')
    if not os.path.isdir(CLONE_DIR):
        os.makedirs(CLONE_DIR)
    copy_tree(htmldir, CLONE_DIR)

    if args.post:
        #-------------------------------------------------------------
        # Post to gh-pages
        #-------------------------------------------------------------

        # cd to directory above where html content was pushed
        os.chdir(CLONE_DIR)
        res, stdout, stderr = get_command_output('touch .nojekyll')
        res1, stdout, stderr1 = get_command_output('git add --all')
        res2, stdout, stderr2 = get_command_output(
            'git commit -am"Pushing version %s to GitHub pages"' % verstr)
        res3, stdout, stderr3 = get_command_output(
            'git push -u origin +gh-pages')
        if res1 + res2 + res3 < 3:
            stderr = stderr1 + stderr2 + stderr3
            print('Something bad happened when attempting to add, commit, '
                  'or push gh-pages content to GitHub - error "%s". Exiting.'
                  % stderr)
            sys.exit(1)
        print('You can inspect the ShakeMap manual and API docs by looking '
              'here: http://usgs.github.io/shakemap/index.html')
    else:
        if not args.clean:
            indexpage = os.path.join(CLONE_DIR, 'index.html')
            print('ShakeMap documentation index: %s' % indexpage)


if __name__ == '__main__':
    desc = 'Create API documentation for ShakeMap and optionally post '\
           'HTML output to GitHub.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-p', '--post', action='store_true', default=False,
                        help='Post documentation to the web.')
    parser.add_argument('-c', '--clean', action='store_true', default=False,
                        help='Clean up local directories containing generated '
                             'HTML content.')
    parser.add_argument('-n', '--nopdf', action='store_true', default=False,
                        help='Suppress PDF building. ')

    pargs = parser.parse_args()
    main(pargs)
