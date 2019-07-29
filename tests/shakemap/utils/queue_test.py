#!/usr/bin/env python

import io
import os
import os.path
import socket
import json
from datetime import datetime, timedelta, timezone
import logging
import time
import shlex
import shutil

import pytest
from shapely.geometry import Polygon

import shakemap.utils.queue as queue
from shakemap.utils.queue import Queue
from shakemap.utils.queue import EventQueue
from shakemap.utils.config import get_config_paths
from shakemap.utils.amps import AmplitudeHandler
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
    assert len(set(config['repeats'][0.0]) ^ set([3600, 7200])) == 0
    assert len(set(config['repeats'][5.0]) ^ set([3600, 7200, 10800])) == 0
    assert config['port'] == 8796
    assert len(set(config['network_delays'].keys()) ^
               set(['ci', 'nc', 'pn'])) == 0
    assert config['network_delays']['ci'] == 300
    assert config['network_delays']['nc'] == 360
    assert config['network_delays']['pn'] == 480


def test_event_queue():
    install_path, _ = get_config_paths()

    db_file = os.path.join(install_path, 'data', 'event_queue.db')
    if os.path.isfile(db_file):
        os.remove(db_file)

    eq = EventQueue(install_path)

    events = eq.getQueuedEvents()
    assert len(events) == 0
    eq.queueEvent('firstevent', ['This', 'is', 'the', 'first', 'event'], 6)
    eq.queueEvent('secondevent', ['This', 'is', 'the', 'second', 'event'], 4.5)
    eq.queueEvent('thirdevent', ['This', 'is', 'the', 'third', 'event'], 6.6)
    events = eq.getQueuedEvents()
    assert events[0][0] == 'thirdevent'
    assert events[0][1] == ['This', 'is', 'the', 'third', 'event']
    assert events[1][0] == 'firstevent'
    assert events[1][1] == ['This', 'is', 'the', 'first', 'event']
    assert events[2][0] == 'secondevent'
    assert events[2][1] == ['This', 'is', 'the', 'second', 'event']
    eq.dequeueEvent('firstevent')
    events = eq.getQueuedEvents()
    assert events[0][0] == 'thirdevent'
    assert events[0][1] == ['This', 'is', 'the', 'third', 'event']
    assert events[1][0] == 'secondevent'
    assert events[1][1] == ['This', 'is', 'the', 'second', 'event']
    eq.dequeueEvent('thirdevent')
    eq.dequeueEvent('secondevent')
    events = eq.getQueuedEvents()
    assert len(events) == 0

    events = eq.getRunningEvents()
    assert len(events) == 0
    eq.insertRunningEvent('firstevent', ['This', 'is', 'the', 'first',
                                         'event'])
    eq.insertRunningEvent('secondevent', ['This', 'is', 'the', 'second',
                                          'event'])
    events = eq.getRunningEvents()
    assert len(events) == 2
    assert events.index(('firstevent', ['This', 'is', 'the', 'first',
                                        'event'])) >= 0
    assert events.index(('secondevent', ['This', 'is', 'the', 'second',
                                         'event'])) >= 0
    eq.deleteRunningEvent('firstevent')
    events = eq.getRunningEvents()
    assert len(events) == 1
    eq.deleteRunningEvent('secondevent')
    events = eq.getRunningEvents()
    assert len(events) == 0

    db_file = eq.db_file
    del eq
    os.remove(db_file)


def test_queue():
    install_path, _ = get_config_paths()

    db_file = os.path.join(install_path, 'data', 'event_queue.db')
    if os.path.isfile(db_file):
        os.remove(db_file)

    pargs = type('Args', (), {})()
    pargs.attached = False
    pargs.break_lock = True
    qq = Queue(pargs)
    qq.logger = qq.getLogger()
    qq.eventQueue = EventQueue(qq.install_path)
    qq.ampHandler = AmplitudeHandler(qq.install_path, qq.data_path)

    assert qq.magnitudeTooSmall({'mag': 2.0,
                                 'lon': -118.25,
                                 'lat': 34.0}) is True
    assert qq.magnitudeTooSmall({'mag': 3.6,
                                 'lon': -118.25,
                                 'lat': 34.0}) is False
    assert qq.magnitudeTooSmall({'mag': 2.0,
                                 'lon': -119.25,
                                 'lat': 35.0}) is True
    assert qq.magnitudeTooSmall({'mag': 3.9,
                                 'lon': -119.25,
                                 'lat': 34.6}) is False
    assert qq.magnitudeTooSmall({'mag': 2.0,
                                 'lon': -129.25,
                                 'lat': 39.0}) is True
    assert qq.magnitudeTooSmall({'mag': 4.1,
                                 'lon': -129.25,
                                 'lat': 39.0}) is False

    dt = datetime.utcnow()

    delta = timedelta(days=180)
    dt_past = dt - delta
    event = {'time': dt_past.strftime(constants.TIMEFMT)}
    assert qq.eventTooOldOrInFuture(event) is False

    delta = timedelta(days=400)
    dt_past = dt - delta
    event = {'time': dt_past.strftime(constants.TIMEFMT)}
    assert qq.eventTooOldOrInFuture(event) is True

    delta = timedelta(minutes=4)
    dt_future = dt + delta
    event = {'time': dt_future.strftime(constants.TIMEFMT)}
    assert qq.eventTooOldOrInFuture(event) is False

    delta = timedelta(minutes=10)
    dt_future = dt + delta
    event = {'time': dt_future.strftime(constants.TIMEFMT)}
    assert qq.eventTooOldOrInFuture(event) is True

    qq.logger.info('Testing the logger')
    fd = open(os.path.join(qq.logpath, 'queue.log'), 'r')
    lines = fd.readlines()
    fd.close()
    assert 'Testing the logger' in lines[-1]

    # Test various paths for origins

    qq.shake_cmds = shlex.split('associate dyfi select assemble -c "Autorun" '
                                'model mapping')
    dt = datetime.utcnow()
    event = {'id': 'test_event',
             'netid': 'xx',
             'network': 'None',
             'lon': 0,
             'lat': 0,
             'depth': 0,
             'mag': 6.6,
             'time': dt.strftime(constants.TIMEFMT),
             'locstring': 'Nowhere',
             'mech': 'ALL',
             'reference': 'Test',
             'productcode': 'test_event'}

    # Test 'test'
    qq.ampHandler.deleteEvent('test_event')
    qq.processOrigin(event, 'test')
    events = qq.eventQueue.getQueuedEvents()
    assert events[0][0] == 'test_event'
    assert events[0][1][0] == 'echo'
    qq.runQueuedEvents()
    events = qq.eventQueue.getRunningEvents()
    assert events[0][0] == 'test_event'
    assert events[0][1][0] == 'echo'
    time.sleep(1)
    deleted = qq.reapChildren()
    assert len(deleted) == 1
    assert deleted[0] == 'test_event'

    qq.config['shake_path'] = 'echo'

    # Test 'cancel'
    qq.ampHandler.deleteEvent('test_event')
    qq.dispatchEvent(event, 'cancel')
    events = qq.eventQueue.getRunningEvents()
    assert events[0][0] == 'test_event'
    assert events[0][1][0] == 'echo'
    time.sleep(1)
    deleted = qq.reapChildren()
    assert len(deleted) == 1
    assert deleted[0] == 'test_event'

    # Test 'shake'

    qq.ampHandler.deleteEvent('test_event')
    qq.processOrigin(event, 'Event added')
    events = qq.eventQueue.getQueuedEvents()
    assert events[0][0] == 'test_event'
    assert events[0][1][0] == 'echo'
    qq.runQueuedEvents()
    events = qq.eventQueue.getRunningEvents()
    assert events[0][0] == 'test_event'
    assert events[0][1][0] == 'echo'
    time.sleep(1)
    deleted = qq.reapChildren()
    assert len(deleted) == 1
    assert deleted[0] == 'test_event'
    shutil.rmtree(os.path.join(qq.data_path, 'test_event'))
    qq.ampHandler.deleteEvent('test_event')

    # Test 'shake' with a network delay

    qq.ampHandler.deleteEvent('test_event')
    event['netid'] = 'ci'
    qq.processOrigin(event, 'Event added')
    events = qq.eventQueue.getQueuedEvents()
    assert len(events) == 0
    test_event = qq.ampHandler.getEvent('test_event')
    assert test_event['repeats'][0] == \
        int(dt.replace(tzinfo=timezone.utc).timestamp()) + 300
    qq.ampHandler.deleteEvent('test_event')

    if os.path.isfile(db_file):
        os.remove(db_file)

    return


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_send_queue()
    test_str_to_seconds()
    test_get_config()
    test_event_queue()
    test_queue()
