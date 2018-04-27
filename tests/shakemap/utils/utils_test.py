#!/usr/bin/env python

import os.path
import tempfile
import textwrap
import sys
import pytest
import copy
import logging

from configobj import ConfigObj

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))

from validate import ValidateError

from shakemap.utils.utils import get_network_name


def test_network():
    assert get_network_name('foo') == 'unknown'
    assert get_network_name(
        'us') == 'USGS National Earthquake Information Center, PDE'


if __name__ == '__main__':
    test_network()
