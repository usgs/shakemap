#!/bin/bash
echo $PATH

VENV=shake
PYVER=3.5


DEPARRAY=(numpy scipy matplotlib jupyter rasterio fiona xlrd xlwt pandas pytables cartopy shapely h5py gdal descartes sphinx configobj pyproj pytest pytest-cov pytest-mpl psutil lxml flake8 pep8-naming openpyxl)

# turn off whatever other virtual environment user might be in
source deactivate

#remove any previous virtual environments called pager
CWD=`pwd`
cd $HOME;
conda remove --name $VENV --all -y
cd $CWD

conda create --name $VENV --yes --channel conda-forge python=$PYVER ${DEPARRAY[*]} -y

# activate the new environment
source activate shake

# do pip installs of those things that are not available via conda.

#grab the bleeding edge for GEM hazardlib.  They have actual releases
#we can resort to if this becomes a problem.
curl --max-time 60 --retry 3 -L https://github.com/gem/oq-engine/archive/master.zip -o openquake.zip
pip -v install --no-deps openquake.zip
rm openquake.zip

#get MapIO
curl --max-time 60 --retry 3 -L https://github.com/usgs/MapIO/archive/master.zip -o mapio.zip
pip -v install --no-deps mapio.zip
rm mapio.zip

#get impactutils
curl --max-time 60 --retry 3 -L https://github.com/usgs/earthquake-impact-utils/archive/master.zip -o impact.zip
pip -v install --no-deps impact.zip
rm impact.zip

#get shakelib
curl --max-time 60 --retry 3 -L https://github.com/usgs/shakelib/archive/master.zip -o shakelib.zip
pip -v install --no-deps shakelib.zip
rm shakelib.zip

pip install sphinx_rtd_theme

# tell the user they have to activate this environment
echo "Type 'source activate shake' to use this new virtual environment."
