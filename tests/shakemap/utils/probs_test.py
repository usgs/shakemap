#!/usr/bin/env python

import os.path

from shakemap.utils.config import get_config_paths
from shakemap.utils.probs import get_probs
from shakemap.coremods.select import validate_config
from configobj import ConfigObj


def test_probs():
    install_path, config_path = get_config_paths()
    config = ConfigObj(os.path.join(install_path, 'config', 'select.conf'))
    validate_config(config, install_path)
    eid = 'us2000cmy3'


if __name__ == '__main__':
    test_probs()
