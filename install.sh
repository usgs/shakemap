#!/bin/bash

unamestr=`uname`
if [ "$unamestr" == 'Linux' ]; then
    prof=~/.bashrc
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    matplotlibdir=~/.config/matplotlib
elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
    prof=~/.bash_profile
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    matplotlibdir=~/.matplotlib
else
    echo "Unsupported environment. Exiting."
    exit
fi

CC_PKG=c-compiler

source $prof

# Name of virtual environment
VENV=shakemap

developer=0
openquake_deps=0
py_ver=3.8
while getopts p:d:q FLAG; do
  case $FLAG in
    p)
        py_ver=$OPTARG
      ;;
    d)
        echo "Installing developer packages."
        developer=1
      ;;
    q)
        echo "Installing full OpenQuake dependencies."
        openquake_deps=1
      ;;
  esac
done

#if [ $py_ver == '3.8' ] && [ "$unamestr" == 'Linux' ]; then
#    echo "WARNING: ShakeMap on Python v3.8 on some versions of Linux "
#    echo "may fail in unexpected ways. We are enforcing the use "
#    echo "of Python v3.7 until this warning disappears."
#    echo ""
#    py_ver=3.7
#fi

echo "Using python version $py_ver"
echo ""

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

    curl -L $mini_conda_url -o miniconda.sh;

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

    # remove the shell script
    rm miniconda.sh
else
    echo "conda detected, installing $VENV environment..."
fi

echo "Installing mamba from conda-forge"

conda install mamba -y -n base -c conda-forge

echo "Installing packages from conda-forge"

# Choose an environment file based on platform
# only add this line if it does not already exist
grep "/etc/profile.d/conda.sh" $prof
if [ $? -ne 0 ]; then
    echo ". $_CONDA_ROOT/etc/profile.d/conda.sh" >> $prof
fi


# Start in conda base environment
echo "Activate base virtual environment"
eval "$(conda shell.bash hook)"                                                
conda activate base

# Remove existing shakemap environment if it exists
conda remove -y -n $VENV --all
conda clean -y --all

# Extra packages to install with dev option
dev_list=(
    "autopep8>=1.5.7"
    "flake8>=3.9.2"
    "pyflakes>=2.3.1"
    "rope>=0.19.0"
    "yapf>=0.31.0"
    "sphinx>=4.0.2"
    "sphinx-argparse>=0.2.5"
)

openquake_list=(
      "decorator>=4.3"
      "django>=3.2"
      "requests>=2.20"
      "setuptools"
      "toml"
)

# Required package list:
package_list=(
      "python=$py_ver"
      "$CC_PKG"
      "cartopy>=0.20"
      "configobj>=5.0.6"
      "cython>=0.29.23"
      "defusedxml>=0.7.1"
      "descartes>=1.1.0"
      "docutils>=0.14"
      "fiona>=1.8.13"
      "gdal>=3.0.2"
      "h5py>=2.10.0"
      "impactutils>=0.8.28"
      "ipython>=7.22.0"
      "libcomcat>=2.0.12"
      "lockfile>=0.12.2"
      "mapio>=0.7.27"
      "matplotlib-base>=3.4.2"
      "numpy==1.20"
      "obspy>=1.2.2"
      "openmp>=8.0.1"
      "pandas>=1.2.5"
      "ps2ff>=1.5.2"
      "psutil>=5.6.7"
      "pyproj>=2.6.1"
      "pytest>=6.2.4"
      "pytest-cov>=2.12.1"
      "python-daemon>=2.3.0"
      "pytest-faulthandler>=2.0.1"
      "pytest-azurepipelines>=0.8.0"
      "pyzmq<20.0"
      "rasterio==1.2.5"
      "scikit-image>=0.16.2"
      "scipy>=1.4.1"
      "shapely>=1.7.1"
      "simplekml>=1.3.5"
      "strec>=2.1.7"
      "versioneer>=0.20"
      "vcrpy>=4.1.1"
)

if [ $developer == 1 ]; then
    package_list=( "${package_list[@]}" "${dev_list[@]}" )
    echo ${package_list[*]}
fi

if [ $openquake_deps == 1 ]; then
    package_list=( "${package_list[@]}" "${openquake_list[@]}" )
    echo ${package_list[*]}
fi

# Create a conda virtual environment
conda config --add channels 'conda-forge'
conda config --add channels 'defaults'
conda config --set channel_priority flexible

echo "Creating the $VENV virtual environment:"
mamba create -y -n $VENV ${package_list[*]}

# Bail out at this point if the conda create command fails.
# Clean up zip files we've downloaded
if [ $? -ne 0 ]; then
    echo "Failed to create conda environment.  Resolve any conflicts, then try again."
    exit 1
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

# The presence of a __pycache__ folder in bin/ can cause the pip
# install to fail... just to be safe, we'll delete it here.
if [ -d bin/__pycache__ ]; then
    rm -rf bin/__pycache__
fi

if [ $developer == 1 ]; then
    pip install sphinx-argparse
fi

# Install OQ from github to get NGA East since it isn't in a release yet.
echo "Installing OpenQuake from github..."
pip install --upgrade --no-dependencies https://github.com/gem/oq-engine/archive/refs/tags/v3.12.0.tar.gz
if [ $? -ne 0 ];then
    echo "Failed to pip install OpenQuake. Exiting."
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
