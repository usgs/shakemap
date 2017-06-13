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
pip -v install https://github.com/gem/oq-hazardlib/archive/master.zip
pip -v install https://github.com/usgs/MapIO/archive/master.zip
pip -v install https://github.com/usgs/earthquake-impact-utils/archive/master.zip
pip install sphinx_rtd_theme

# tell the user they have to activate this environment
echo "Type 'source activate shake' to use this new virtual environment."
