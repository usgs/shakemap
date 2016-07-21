#!/bin/sh

#where should .rst files, Makefile, _build folder be written?
SPHINX_DIR=$HOME/sphinx

#where should the temporary clone of the shakemap repo live?
CLONE_DIR=$HOME/tmpshake

#get the human friendly version of the ShakeMap version
VERSION=`git tag`

#we should already be in the repo folder where the code lives...
REPO_DIR=`pwd`

#what is the package called and who are the authors
PACKAGE="shakemap"
AUTHORS='Bruce Worden, Eric Thompson, Mike Hearne'

# #clone the repository and cd to it
git clone https://github.com/mhearne-usgs/shakemap.git $CLONE_DIR
cd $CLONE_DIR

#clear out everything that's there now
git rm -rf *

#checkout an orphan branch called gh-pages
#git checkout --orphan gh-pages

#Run the sphinx command that creates the initial set of .rst files, Makefile, and conf.py.
sphinx-apidoc -o $SPHINX_DIR -f -l -F -H $PACKAGE -A "$AUTHORS" -V $VERSION $REPO_DIR

#Modify the conf.py file to include path to repository, read-the-docs theme.
echo "sys.path.insert(0, os.path.abspath('${REPO_DIR}'))" >> $SPHINX_DIR/conf.py
echo "import sphinx_rtd_theme" >> $SPHINX_DIR/conf.py
echo "html_theme = 'sphinx_rtd_theme'" >> $SPHINX_DIR/conf.py
echo "html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]" >> $SPHINX_DIR/conf.py

#Go to the sphinx directory and build the html
cd $SPHINX_DIR
make html

#copy the generated html to the gh-pages repository

#mkdir $CLONE_DIR
cp -R $SPHINX_DIR/_build/html/* $CLONE_DIR

#go to the gh-pages repository, add the new files, commit, and upload to GitHub
cd $CLONE_DIR
touch $CLONE_DIR/.nojekyll
git add --all
git commit -am"Pushing version ${VERSION} to GitHub pages"
git push -u origin +gh-pages

#go back to where we started and then clean up.
rm -rf $SPHINX_DIR
