#!/bin/bash
echo $PATH

DEPS="numpy scipy matplotlib jupyter rasterio fiona xlrd xlwt pandas pytables basemap basemap-data-hires shapely h5py basemap-data-hires gdal==1.11.4 descartes paramiko sphinx configobj pyproj pytest pytest-cov pytest-mpl psutil lxml flake8 pep8-naming"

if [ "$#" -le 1 ]; then
    # turn off whatever other virtual environment user might be in
    source deactivate
    
    # remove any previous virtual environments called shake
    conda remove --name shake --all -y
    
    # create a new virtual environment called shake with the below list of dependencies installed into it
    conda create --name shake --yes --channel conda-forge python=3.5 $DEPS -y
else
    conda install --yes --channel conda-forge python=3.5 $DEPS -y
fi

# activate the new environment
source activate shake

# do pip installs of those things that are not available via conda.
pip -v install git+git://github.com/gem/oq-hazardlib.git
pip install git+git://github.com/usgs/MapIO.git
pip install git+git://github.com/usgs/earthquake-impact-utils.git
pip install sphinx_rtd_theme

# tell the user they have to activate this environment
echo "Type 'source activate shake' to use this new virtual environment."
