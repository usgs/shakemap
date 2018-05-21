#!/bin/bash

unamestr=`uname`
if [ "$unamestr" == 'Linux' ]; then
    source ~/.bashrc
elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
    source ~/.bash_profile
fi

echo "Path:"
echo $PATH

VENV=shakemap

# Is the reset flag set?
reset=0
while getopts r FLAG; do
  case $FLAG in
    r)
        reset=1
        
      ;;
  esac
done

# Get uname to determine OS
unamestr=`uname`


# create a matplotlibrc file with the non-interactive backend "Agg" in it.
if [ ! -d ~/.config/matplotlib ]; then
    mkdir -p ~/.config/matplotlib
    echo "backend : Agg" > ~/.config/matplotlib/matplotlibrc
    echo "NOTE: A non-interactive matplotlib backend (Agg) has been set for this user."
fi

# Is conda installed?
conda --version
if [ $? -ne 0 ]; then
    echo "No conda detected, installing miniconda..."
    if [ "$unamestr" == 'Linux' ]; then
	mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
	mini_conda_url=https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
    else
	echo "Unsupported environment. Exiting."
	exit
    fi

    curl $mini_conda_url -o miniconda.sh;
    echo "Install directory: $HOME/miniconda"

    bash miniconda.sh -f -b -p $HOME/miniconda

    # Need this to get conda into path
    . $HOME/miniconda/etc/profile.d/conda.sh
else
    echo "conda detected, installing $VENV environment..."
fi

echo "PATH:"
echo $PATH
echo ""


# Choose an environment file based on platform
if [ "$unamestr" == 'Linux' ]; then
    env_file=environment_linux.yml
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
    env_file=environment_osx.yml
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bash_profile
else
    echo "Unsupported environment. Exiting."
    exit
fi

# If the user has specified the -r (reset) flag, then create an
# environment based on only the named dependencies, without
# any versions of packages specified.
if [ $reset == 1 ]; then
    echo "Ignoring platform, letting conda sort out dependencies..."
    env_file=environment.yml
fi

# Start in conda base environment
echo "Activate base virtual environment"
conda activate base

# Create a conda virtual environment
echo "Creating the $VENV virtual environment:"
conda env create -f $env_file --force

# Bail out at this point if the conda create command fails.
# Clean up zip files we've downloaded
if [ $? -ne 0 ]; then
    echo "Failed to create conda environment.  Resolve any conflicts, then try again."
    exit
fi


# Activate the new environment
echo "Activating the $VENV virtual environment"
conda activate $VENV

# This package
echo "Installing shakemap..."
pip install -e .

# Install default profile
#python bin/sm_profile -c default -a

# Tell the user they have to activate this environment
echo "Type 'conda activate $VENV' to use this new virtual environment."
