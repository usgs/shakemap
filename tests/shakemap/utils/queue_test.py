#!/usr/bin/env python

import io
import os
import os.path
import socket
import json
from datetime import datetime, timedelta
import logging
import subprocess
import time

import pytest
from shapely.geometry import Polygon

import shakemap.utils.queue as queue
from shakemap.utils.config import get_config_paths
from shakelib.rupture import constants


def get_dummy_logger(name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    log_capture_string = io.StringIO()
    handler = logging.StreamHandler(log_capture_string)
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger, log_capture_string


def test_send_queue():
    install_path, _ = get_config_paths()
    config = queue.get_config(install_path)
    rsocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    rsocket.bind(('', config['port']))
    rsocket.listen(1)

    queue.send_queue('test', {'id': 'testid'}, config['port'])

    (clientsocket, address) = rsocket.accept()
    data = clientsocket.recv(queue.MAX_SIZE)
    clientsocket.close()
    cmd = json.loads(data.decode('utf-8'))
    assert cmd['type'] == 'test'
    assert cmd['data']['id'] == 'testid'

    # Send a too-large data segment
    bigdata = [x for x in range(queue.MAX_SIZE)]
    with pytest.raises(RuntimeError):
        queue.send_queue('test', bigdata)


def test_str_to_seconds():
    assert queue.str_to_seconds('1m') == 60
    assert queue.str_to_seconds('1h') == 3600
    assert queue.str_to_seconds('1d') == 86400
    assert queue.str_to_seconds('1y') == 31536000
    assert queue.str_to_seconds('1') == 60
    assert queue.str_to_seconds('-1') == -1
    with pytest.raises(ValueError):
        queue.str_to_seconds('kkkaw')


def test_get_config():
    install_path, _ = get_config_paths()
    config = queue.get_config(install_path)
    assert len(set(config['servers']) ^ set(['localhost'])) == 0
    assert config['max_process_time'] == 600
    assert config['old_event_age'] == queue.str_to_seconds('1y')
    assert config['future_event_age'] == queue.str_to_seconds('5m')
    assert config['minmag'] == 4.0
    p1 = Polygon([(-116.75, 33.50), (-118.25, 33.50), (-120.25, 34.33),
                  (-120.25, 34.75), (-116.75, 34.75), (-116.75, 33.50)])
    assert config['boxes']['01_my_box']['mag'] == 3.5
    assert config['boxes']['01_my_box']['poly'] == p1
    p2 = Polygon([(-117.75, 34.50), (-119.25, 34.50), (-121.25, 34.33),
                  (-121.25, 34.75), (-117.75, 34.75), (-117.75, 34.50)])
    assert config['boxes']['02_my_box2']['mag'] == 3.8
    assert config['boxes']['02_my_box2']['poly'] == p2
    assert len(set(config['repeats'][0.0]) ^ set([60, 120])) == 0
    assert len(set(config['repeats'][5.0]) ^ set([60, 120, 180])) == 0
    assert config['port'] == 8796


def test_magnitude_too_small():
    install_path, _ = get_config_paths()
    config = queue.get_config(install_path)
    assert queue.magnitude_too_small(2.0, -118.25, 34.0, config) is True
    assert queue.magnitude_too_small(3.6, -118.25, 34.0, config) is False
    assert queue.magnitude_too_small(2.0, -119.25, 35.0, config) is True
    assert queue.magnitude_too_small(3.9, -119.25, 34.6, config) is False
    assert queue.magnitude_too_small(2.0, -129.25, 39.0, config) is True
    assert queue.magnitude_too_small(4.1, -129.25, 39.0, config) is False


def test_event_too_old_or_in_future():
    install_path, _ = get_config_paths()
    config = queue.get_config(install_path)
    dt = datetime.utcnow()

    delta = timedelta(days=180)
    dt_past = dt - delta
    event = {'time': dt_past.strftime(constants.TIMEFMT)}
    assert queue.event_too_old_or_in_future(event, config) is False

    delta = timedelta(days=400)
    dt_past = dt - delta
    event = {'time': dt_past.strftime(constants.TIMEFMT)}
    assert queue.event_too_old_or_in_future(event, config) is True

    delta = timedelta(minutes=4)
    dt_future = dt + delta
    event = {'time': dt_future.strftime(constants.TIMEFMT)}
    assert queue.event_too_old_or_in_future(event, config) is False

    delta = timedelta(minutes=10)
    dt_future = dt + delta
    event = {'time': dt_future.strftime(constants.TIMEFMT)}
    assert queue.event_too_old_or_in_future(event, config) is True


def test_get_logger():
    install_path, _ = get_config_paths()
    logpath = os.path.join(install_path, 'logs')
    logger = queue.get_logger(logpath, False)
    logger.info('Testing the logger')
    fd = open(os.path.join(logpath, 'queue.log'), 'r')
    lines = fd.readlines()
    fd.close()
    assert 'Testing the logger' in lines[-1]


def test_dispatch_event():
    # Test 'test'
    logger, logstring = get_dummy_logger('test_dispatch')
    children = {}
    config = {}
    shake_config = {'autorun_modules': 'associate dyfi select assemble '
                                       '-c "Autorun" model mapping'}
    queue.dispatch_event('Testing dispatch', logger, children, 'test', config,
                         shake_config)
    returncode = children['Testing dispatch']['popen'].wait()
    assert returncode == 0
    assert 'Testing event Testing dispatch' in logstring.getvalue()

    # Test 'cancel'
    children = {}
    config = {'cancel_command': 'echo "Testing cancel"',
              'shake_path': '/dev/null'}
    queue.dispatch_event('Testing dispatch', logger, children, 'cancel',
                         config, shake_config)
    returncode = children['Testing dispatch']['popen'].wait()
    assert returncode == 0
    assert 'Canceling event Testing dispatch' in logstring.getvalue()

    # Test 'shake'
    children = {}
    config = {'shake_command': 'echo "Testing shake"',
              'shake_path': '/dev/null'}
    queue.dispatch_event('Testing dispatch', logger, children, 'Event added',
                         config, shake_config)
    returncode = children['Testing dispatch']['popen'].wait()
    assert returncode == 0
    assert 'Running event Testing dispatch due to action "Event added"' in \
        logstring.getvalue()
    logstring.close()
    return


def test_reap_children():
    logger, logstring = get_dummy_logger('test_reap_children')
    config = {'max_process_time': 10}
    p = subprocess.Popen(['echo', 'Testing reap'])
    children = {}
    children['Testing reap'] = {'popen': p, 'start_time': time.time()}
    time.sleep(0.5)
    queue.reap_children(children, config, logger)
    assert len(list(children.keys())) == 0
    assert 'Reaped child for event Testing reap' in logstring.getvalue()

    config = {'max_process_time': 1}
    p = subprocess.Popen(['sleep', '10'])
    children['Testing kill'] = {'popen': p, 'start_time': time.time()}
    time.sleep(2.0)
    queue.reap_children(children, config, logger)
    assert len(list(children.keys())) == 0
    assert 'Event Testing kill taking too long, killing' in \
        logstring.getvalue()
    logstring.close()

    return


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_send_queue()
    test_str_to_seconds()
    test_get_config()
    test_magnitude_too_small()
    test_event_too_old_or_in_future()
    test_get_logger()
    test_dispatch_event()
    test_reap_children()
