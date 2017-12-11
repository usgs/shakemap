#!/bin/bash

echo "Path:"
echo $PATH

VENV=shakemap

# Is conda installed?
conda=$(which conda)
if [ ! "$conda" ] ; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O miniconda.sh;
    bash miniconda.sh -f -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
fi

unamestr=`uname`
if [ "$unamestr" == 'Linux' ]; then
    env_file=environment_linux.yml
elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
    env_file=environment_osx.yml
fi

# Turn off whatever other virtual environment user might be in
source deactivate

# Download dependencies not in conda or pypi
curl --max-time 60 --retry 3 -L \
    https://github.com/usgs/earthquake-impact-utils/archive/master.zip -o impact-utils.zip
curl --max-time 60 --retry 3 -L \
    https://github.com/usgs/libcomcat/archive/master.zip -o libcomcat.zip
curl --max-time 60 --retry 3 -L \
    https://github.com/usgs/MapIO/archive/master.zip -o mapio.zip


# Create a conda virtual environment
echo "Creating the $VENV virtual environment:"
conda env create -f $env_file --force

if [ $? -ne 0 ]; then
    echo "Failed to create conda environment.  Resolve any conflicts, then try again."
    exit
fi


# Activate the new environment
echo "Activating the $VENV virtual environment"
source activate $VENV

# Install OpenQuake -- note that I have pulled this out of environment.yml
# because the requirements are too narrow to work with our other dependencies,
# but the openquake.hazardlib tests pass with this environment. We need to
# remember to check this when we change the environemnt.yml file though.
conda install -y --no-deps -c conda-forge openquake.engine


# Clean up downloaded packages
rm impact-utils.zip
rm libcomcat.zip
rm mapio.zip


# This package
echo "Installing shakemap..."
pip install -e .

# Install default profile
#python bin/sm_profile -c default -a

# Tell the user they have to activate this environment
echo "Type 'source activate $VENV' to use this new virtual environment."
