import os

#
# This is needed here so that the matplotlib backend gets
# set before any other imports of matplotlib
#
import matplotlib
matplotlib.use('Agg')


def pytest_configure(config):
    #
    # This tells get_config_paths() (shakemap.utils.config) to
    # return paths into the testing part of the repo
    #
    os.environ['CALLED_FROM_PYTEST'] = 'True'


def pytest_unconfigure(config):
    del os.environ['CALLED_FROM_PYTEST']
