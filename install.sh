#!/bin/bash

unamestr=`uname`
if [ "$unamestr" == 'Linux' ]; then
    prof=~/.bashrc
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    matplotlibdir=~/.config/matplotlib
    channel_url=ftp://ftpext.usgs.gov/pub/cr/co/golden/emthompson/shakemap-linux.tar
    channel=shakemap-linux
elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
    prof=~/.bash_profile
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    matplotlibdir=~/.matplotlib
    channel_url=ftp://ftpext.usgs.gov/pub/cr/co/golden/emthompson/shakemap-osx.tar
    channel=shakemap-osx
else
    echo "Unsupported environment. Exiting."
    exit
fi

source $prof

# Name of virtual environment
VENV=shakemap

# Is the reset flag set? If not, use the frozen conda channel to install.
# Otherwise, try to install from conda forge.
reset=0
while getopts r FLAG; do
  case $FLAG in
    r)
        reset=1
        
      ;;
  esac
done


# create a matplotlibrc file with the non-interactive backend "Agg" in it.
if [ ! -d "$matplotlibdir" ]; then
    mkdir -p $matplotlibdir
    # if mkdir fails, bow out gracefully
    if [ $? -ne 0 ];then
        echo "Failed to create matplotlib configuration file. Exiting."
        exit 1
    fi
fi
matplotlibrc=$matplotlibdir/matplotlibrc
if [ ! -e "$matplotlibrc" ]; then
    echo "backend : Agg" > "$matplotlibrc"
    echo "NOTE: A non-interactive matplotlib backend (Agg) has been set for this user."
elif grep -Fxq "backend : Agg" $matplotlibrc ; then
    :
elif [ ! grep -Fxq "backend" $matplotlibrc ]; then
    echo "backend : Agg" >> $matplotlibrc
    echo "NOTE: A non-interactive matplotlib backend (Agg) has been set for this user."
else
    sed -i '' 's/backend.*/backend : Agg/' $matplotlibrc
    echo "###############"
    echo "NOTE: $matplotlibrc has been changed to set 'backend : Agg'"
    echo "###############"
fi


# Is conda installed?
conda --version
if [ $? -ne 0 ]; then
    echo "No conda detected, installing miniconda..."

    command -v curl >/dev/null 2>&1 || { echo >&2 "Script requires curl but it's not installed. Aborting."; exit 1; }

    curl $mini_conda_url -o miniconda.sh;

    # if curl fails, bow out gracefully
    if [ $? -ne 0 ];then
        echo "Failed to download miniconda installer shell script. Exiting."
        exit 1
    fi
    
    echo "Install directory: $HOME/miniconda"

    bash miniconda.sh -f -b -p $HOME/miniconda

    # if miniconda.sh fails, bow out gracefully
    if [ $? -ne 0 ];then
        echo "Failed to run miniconda installer shell script. Exiting."
        exit 1
    fi
    
    . $HOME/miniconda/etc/profile.d/conda.sh
else
    echo "conda detected, installing $VENV environment..."
fi

if [ $reset == 0 ]; then
    # Download frozen channel
    echo "Downloading shakemap channel..."
    curl $channel_url -o $channel.tar
    # if curl fails, bow out gracefully
    if [ $? -ne 0 ];then
	echo "Failed to download channel. Exiting."
	exit 1
    fi

    # Un tar the channel
    tar -xvf $channel.tar
    if [ $? -ne 0 ];then
	echo "Failed to extract channel. Exiting."
	exit 1
    fi
else
    echo "Installing packages from conda-forge"
fi

# Choose an environment file based on platform
# only add this line if it does not already exist
grep "/etc/profile.d/conda.sh" $prof
if [ $? -ne 0 ]; then
    echo ". $_CONDA_ROOT/etc/profile.d/conda.sh" >> $prof
fi

env_file=environment.yml


# Start in conda base environment
echo "Activate base virtual environment"
conda activate base

# Remove existing shakemap environment if it exists
conda remove -y -n $VENV --all


# Package list:
package_list='
      python=3.5
      amptools
      basemap
      cartopy
      cython
      defusedxml
      descartes
      docutils
      configobj
      fiona
      gcc
      gdal
      h5py
      impactutils
      libcomcat
      lockfile
      mapio
      matplotlib<=2.3
      numexpr
      numpy
      obspy
      openquake.engine
      pandas
      ps2ff
      psutil
      pyproj
      pytest
      pytest-cov
      python-daemon
      pytest-faulthandler
      scikit-image
      scipy
      shapely
      simplekml
      strec
      versioneer 
      vcrpy
'

# Create a conda virtual environment
echo "Creating the $VENV virtual environment:"
if [ $reset == 0 ]; then
    conda create -y --override-channels -n $VENV \
          -c file://$PWD/$channel $package_list
else
    conda create -y -n $VENV -c conda-forge $package_list
fi


# Bail out at this point if the conda create command fails.
# Clean up zip files we've downloaded
if [ $? -ne 0 ]; then
    echo "Failed to create conda environment.  Resolve any conflicts, then try again."
    exit
fi


# Activate the new environment
echo "Activating the $VENV virtual environment"
conda activate $VENV

# if conda activate fails, bow out gracefully
if [ $? -ne 0 ];then
    echo "Failed to activate ${VENV} conda environment. Exiting."
    exit 1
fi

# upgrade pip, mostly so pip doesn't complain about not being new...
pip install --upgrade pip

# if pip upgrade fails, complain but try to keep going
if [ $? -ne 0 ];then
    echo "Failed to upgrade pip, trying to continue..."
    exit 1
fi

# This package
echo "Installing ${VENV}..."
pip install --no-deps -e .

# if pip install fails, bow out gracefully
if [ $? -ne 0 ];then
    echo "Failed to pip install this package. Exiting."
    exit 1
fi

# Tell the user they have to activate this environment
echo "Type 'conda activate $VENV' to use this new virtual environment."
