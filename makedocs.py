#!/usr/bin/env python

import argparse
import subprocess
import os.path
import shutil
import sys
import re

from shakemap.utils.misc import getCommandOutput

DEFAULT_TAG = '0.1'


def main(args):

    # where should .rst files, Makefile, _build folder be written?
    SPHINX_DIR = os.path.join(os.path.expanduser('~'), '__api-doc')
    print("Sphinx dir: %s" %SPHINX_DIR)
    shutil.rmtree(SPHINX_DIR, ignore_errors=True)

    # where should the temporary clone of the shakemap gh-pages repo live?
    TOP_DIR=os.path.join(os.path.expanduser('~'),'__shake-doc')
    CLONE_DIR=os.path.join(TOP_DIR,'html')
    print("Clone dir: %s" %CLONE_DIR)
    shutil.rmtree(TOP_DIR, ignore_errors=True)
    # where is the repository where this script is running?
    REPO_DIR = os.path.dirname(os.path.abspath(__file__))  # where is this script?
    print("Repo dir: %s" %REPO_DIR)

    # get the human friendly version of the ShakeMap version
    res, verstr, stderr = getCommandOutput('git describe --tags')
    verstr = verstr.decode().strip()
    if not len(verstr):
        verstr = DEFAULT_TAG
    else:
        if verstr[0] == 'v':
            verstr = verstr[1:]
        spl = re.findall(r"[\w']+", verstr)
        verstr = "%s.%s.%s" %(spl[0], spl[1], spl[2])

    # what is the package called and who are the authors
    PACKAGE = "shakemap"
    AUTHORS = 'Bruce Worden, Eric Thompson, Mike Hearne'

    # find the make command on this system
    res, stdout, stderr = getCommandOutput('which make')
    if not res:
        print('Could not find the "make" command on your system. Exiting.')
        sys.exit(1)
    make_cmd = stdout.decode().strip()

    try:
        # clone the repository
        if args.post:
            sys.stderr.write('Cloning shakemap gh-pages branch...\n')
            if os.path.isdir(CLONE_DIR):
                shutil.rmtree(CLONE_DIR)
            clonecmd = 'git clone -b gh-pages https://github.com/usgs/shakemap.git %s' % CLONE_DIR
            res,stdout,stderr = getCommandOutput(clonecmd)
            if not res:
                raise Exception('Could not clone gh-pages branch.')

        # remove the current apidocs folder from that repo
        apidocfolder = os.path.join(CLONE_DIR, 'apidoc')
        # this should probably not have 'html' in the path
        if os.path.isdir(apidocfolder):
            shutil.rmtree(apidocfolder)

        # run the make command to build the shakemap manual
        sys.stderr.write('Building shakemap manual (HTML)...\n')
        makefile = os.path.join(REPO_DIR, 'doc', 'Makefile')
        os.chdir(os.path.join(REPO_DIR, 'doc'))
        manualcmd = '%s html' % make_cmd
        res,stdout,stderr = getCommandOutput(manualcmd)
        res,stdout,stderr = getCommandOutput(manualcmd)
        if not res:
            raise Exception('Could not build the ShakeMap manual - error "%s".' % stderr)
        
        # run the sphinx api doc command
        sys.stderr.write('Building shakemap API documentation (REST)...\n')
        package_dir = os.path.join(REPO_DIR,'shakemap')
        sphinx_cmd = 'sphinx-apidoc -o %s -f -l -e -d 12 -F -H %s -A "%s" -V %s %s' % (SPHINX_DIR,PACKAGE,AUTHORS,verstr,package_dir)
        print("\n'%s'\n" %sphinx_cmd)
        res,stdout,stderr = getCommandOutput(sphinx_cmd)
        if not res:
            raise Exception('Could not build ShakeMap API documentation - error "%s".' % stderr)

        # run the make command to build the shakemap manual (pdf version)
        if not args.nopdf:
            sys.stderr.write('Building shakemap manual (PDF)...\n')
            os.chdir(os.path.join(REPO_DIR,'doc'))
            manualcmd = '%s latexpdf' % make_cmd
            res,stdout,stderr = getCommandOutput(manualcmd)
            if not res:
                raise Exception('Could not build the PDF version of the ShakeMap manual - error "%s".' % stderr)

        # this has created a conf.py file and a Makefile.
        
        # Edit the conf.py file to include the theme.
        fname = os.path.join(SPHINX_DIR,'conf.py')
        f = open(fname,'at')
        f.write("sys.path.insert(0, os.path.abspath('%s'))\n" % (REPO_DIR))

        f.write("import sphinx_rtd_theme\n")
        f.write("html_theme = 'sphinx_rtd_theme'\n")
        f.write("html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]\n")

        f.write("html_theme_options = {\n")
        f.write("    'collapse_navigation': False,\n")
#        f.write("    'navigation_depth': 12,\n")
        f.write("}\n")

        # This line is needed to inclue __init__ methods in documentation
        f.write("autoclass_content = 'both'\n")
        f.write("html_show_copyright = False\n")
        f.close()

        # Go to the sphinx directory and build the html
        sys.stderr.write('Building shakemap API documentation (HTML)...\n')
        os.chdir(SPHINX_DIR)
        res,stdout,stderr = getCommandOutput('%s html' % make_cmd)
        if not res:
            raise Exception('Could not build HTML for API documentation. - error "%s"' % stderr)
        
        # copy the generated content to the gh-pages branch we created earlier
        htmldir = os.path.join(SPHINX_DIR,'_build','html')
        if os.path.isdir(apidocfolder):
            shutil.rmtree(apidocfolder)
        shutil.copytree(htmldir, apidocfolder)

        if args.post:
            # cd to directory above where html content was pushed
            htmldir = os.path.join(CLONE_DIR) #this should probably not have 'html' in it
            os.chdir(htmldir)
            res,stdout,stderr = getCommandOutput('touch .nojekyll')
            res1,stdout,stderr1 = getCommandOutput('git add --all')
            res2,stdout,stderr2 = getCommandOutput('git commit -am"Pushing version %s to GitHub pages"' % verstr)
            res3,stdout,stderr3 = getCommandOutput('git push -u origin +gh-pages')
            if res1+res2+res3 < 3:
                stderr = stderr1+stderr2+stderr3
                print('Something bad happened when attempting to add, commit, or push gh-pages content to GitHub - error "%s". Exiting.' % stderr)
                sys.exit(1)
            print('You can inspect the ShakeMap manual and API docs by looking here: http://usgs.github.io/shakemap/index.html')
        else:
            if not args.clean:
                indexpage = os.path.join(CLONE_DIR,'index.html') #this should not have 'html' in it
                print('You can inspect the ShakeMap manual and API docs by looking here: %s' % indexpage)
    except Exception as error:
        print('Error detected: %s' % str(error))
    finally:
        if args.clean:
            print('Cleaning up %s and %s created directories...' %
                  (CLONE_DIR, SPHINX_DIR))
            shutil.rmtree(CLONE_DIR)
            shutil.rmtree(SPHINX_DIR)


if __name__ == '__main__':
    desc = '''Create API documentation for ShakeMap and optionally post HTML output to GitHub.
    '''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-p', '--post', action='store_true', default=False,
                        help='Post content to the web')
    parser.add_argument('-c', '--clean', action='store_true', default=False,
                        help='Clean up local directories containing generated HTML content')
    parser.add_argument('-n', '--nopdf', action='store_true', default=False,
                        help='Suppress PDF building. ')

    pargs = parser.parse_args()
    main(pargs)
