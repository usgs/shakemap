import os

def pytest_configure(config):
    os.environ['CALLED_FROM_PYTEST'] = 'True'

def pytest_unconfigure(config):
    del os.environ['CALLED_FROM_PYTEST']
