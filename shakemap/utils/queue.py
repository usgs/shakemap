import os
import os.path
import socket
import json
from collections import OrderedDict
import time as time
from datetime import datetime, timezone
import logging
from logging.handlers import TimedRotatingFileHandler
import subprocess
import shlex

from configobj import ConfigObj
from validate import Validator
from shapely.geometry import (Polygon,
                              Point)

from shakemap.utils.config import (config_error,
                                   get_configspec)
from shakelib.rupture import constants

MAX_SIZE = 4096


def send_queue(command, data, port=9755):
    """
    Send a command and its data to the queue process.

    Args:
        command (str): A valid command for the queue (e.g., 'origin').
        data (JSON serializable): The data associated with the command.
            Could be int, float, string, dictionary, list, etc. Must be
            JSON serializable. Must be less than MAX_SIZE bytes when
            serialized.

    Returns:
        nothing: Nothing

    Raises:
        RuntimeError: If the serialized data is larger than MAX_SIZE.
        OSError: If there is a problem with the socket or connection.
        TypeError: If the supplied data is not JSON serializable.
    """
    qdata = {'type': command, 'data': data}
    qstr = json.dumps(qdata).encode('utf8')
    if len(qstr) > MAX_SIZE:
        raise RuntimeError('Data portion too large')

    csocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    csocket.connect(('127.0.0.1', port))
    csocket.send(qstr)
    csocket.close()
    return


def str_to_seconds(tstring):
    """ Convert time strings to seconds. Strings can be of the
    form:
        <int>   (ninutes)
        <int>m  (minutes)
        <int>h  (hours)
        <int>d  (days)
        <int>y  (years)

    Args:
        tstring (str): An integer followed by an (optional)
                       'm', 'h', 'd', 'y'.

    Returns
        int: The number of seconds represented by the input string. If
        the string is unadorned and represents a negative number, -1
        is returned.

    Raises:
        ValueError: If the value cannot be converted to seconds.
    """
    if tstring.endswith('m'):
        secs = 60 * int(tstring.replace('m', ''))
    elif tstring.endswith('h'):
        secs = 60 * 60 * int(tstring.replace('h', ''))
    elif tstring.endswith('d'):
        secs = 24 * 60 * 60 * int(tstring.replace('d', ''))
    elif tstring.endswith('y'):
        secs = 365 * 24 * 60 * 60 * int(tstring.replace('y', ''))
    else:
        secs = 60 * int(tstring)
        if secs < 0:
            secs = -1

    return secs


def parse_config(config):
    """ Parse the config object to get usable data.

    Args:
        config (ConfigObj object): The result of parsing the file's
                                   configuration file.

    Returns:
        dict: A cleaned up version of the input.
    """
    if '' in config['servers']:
        config['servers'].remove('')
    if 'localhost' not in config['servers']:
        config['servers'].append('localhost')

    config['old_event_age'] = str_to_seconds(config['old_event_age'])
    config['future_event_age'] = str_to_seconds(config['future_event_age'])
    config['associate_interval'] = str_to_seconds(config['associate_interval'])
    config['max_trigger_wait'] = str_to_seconds(config['max_trigger_wait'])

    boxes = OrderedDict()
    for key, value in config['boxes'].items():
        coords = [float(x) for x in value.split(',')]
        mag = coords.pop(0)
        # The config is in lat, lon order, but wee want lon, lat
        # so we do the odd indices first
        coords2 = list(zip(coords[1::2], coords[0::2]))
        coords2.append(coords2[0])
        boxes[key] = {'mag': mag, 'poly': Polygon(coords2)}
    config['boxes'] = boxes

    repeats = {}
    for key, value in config['repeats'].items():
        tlist = [str_to_seconds(x) for x in value]
        repeats[float(key)] = tlist
    config['repeats'] = repeats

    return config


def get_config(install_path):
    """Read the config and get it into a usable state.

    Args:
        install_path (str): The install path of the current profile.

    Returns:
        dict: A dictionary of configuration data.
    """
    config_file = os.path.join(install_path, 'config', 'queue.conf')
    configspec = get_configspec('queue')
    config = ConfigObj(config_file, configspec=configspec)
    results = config.validate(Validator())
    if not isinstance(results, bool) or not results:
        config_error(config, results)
    config = parse_config(config.dict())
    return config


def magnitude_too_small(mag, lon, lat, config):
    """ Return False if the magnitude is greater than the threshold
    magnitude of the first metro box within which it falls, or
    the global minmag if it does not fall within a box; return
    true otherwise.

    Args:
        mag (float): event magnitude
        lon (float): event longitude
        lat (float): event latitude
        config (dict): the configuration structure

    Returns:
        bool: True if the event is too small to process; False otherwise.
    """
    pt = Point((lon, lat))
    for boxname in sorted(config['boxes']):
        boxdict = config['boxes'][boxname]
        if pt.within(boxdict['poly']):
            if mag >= boxdict['mag']:
                return False
            else:
                return True
    #
    # Not in any boxes
    #
    if mag >= config['minmag']:
        return False

    return True


def event_too_old_or_in_future(event, config):
    """ Return True if the event is too old or too far in the future to
    process; return False otherwise.

    Args:
        event (dict): The event data structure.
        config (dict): The configuration data structure

    Returns:
        bool: True if the event is older than old_event_age or is
        more than future_event_age in the future; returns False otherwise.
    """

    current_time = time.time()
    try:
        event_time = datetime.strptime(event['time'], constants.TIMEFMT).\
            replace(tzinfo=timezone.utc).timestamp()
    except ValueError:
        event_time = datetime.strptime(event['time'], constants.ALT_TIMEFMT).\
            replace(tzinfo=timezone.utc).timestamp()
    if config['old_event_age'] >= 0 and \
       event_time + config['old_event_age'] < current_time:
        return True
    if config['future_event_age'] >= 0 and \
       event_time - config['future_event_age'] > current_time:
        return True

    return False


def get_logger(logpath, attached):
    """Set up a logger for this process.

    Args:
        logpath (str): Path to the directory into which to put the logfile.

    Returns:
        logging.logger: An instance of a logger.
    """
    if not os.path.isdir(logpath):
        os.makedirs(logpath)
    logger = logging.getLogger('queue_logger')
    logger.setLevel(logging.INFO)
    if attached:
        handler = logging.StreamHandler()
    else:
        logfile = os.path.join(logpath, 'queue.log')
        handler = TimedRotatingFileHandler(logfile,
                                           when='midnight',
                                           backupCount=30)
    formatter = logging.Formatter(
            fmt='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.propagate = False
    return logger


def dispatch_event(eventid, logger, children, action, config, shake_config):
    """ Start a subprocess to run the specified event and add its data
    to the children structure.

    Args:
        eventid (str): The event ID to process.
        logger (logger): The process logger.
        children (dict): The structure holding info on the child processes.
        action (str): 'cancel', 'test', or some other string. 'cancel' starts
                      the cancel process, 'test' starts the process
                      'echo eventid'.  And any other string starts the
                      shake process. See the configuration file 'queue.conf'
                      for the exact commands that will be run.
        config (dict): The configuration dictionary.
        shake_config (dict): The parsed shake.conf config file.

    Returns:
        nothing: Nothing.
    """
    if action == 'cancel':
        logger.info('Canceling event %s' % eventid)
        cmd = config['cancel_command'].replace('shake', config['shake_path'])
        cmd = cmd.replace('<EVID>', eventid)
        cmd = cmd.split()
        p = subprocess.Popen(cmd)
    elif action == 'test':
        logger.info('Testing event %s' % eventid)
        p = subprocess.Popen(['echo', eventid])
    else:
        logger.info('Running event %s due to action "%s"' % (eventid, action))
        #
        # Add the action as the assemble/augment comment, or replace the
        # comment if it is already there.
        #
        shake_cmds = shlex.split(shake_config['autorun_modules'])
        for ix, shcmd in enumerate(shake_cmds):
            if shcmd not in ['assemble', 'augment']:
                continue
            if len(shake_cmds) == ix + 1:  # This shouldn't happen
                shake_cmds.append('-c')
                shake_cmds.append('"%s"' % action)
            elif shake_cmds[ix + 1] == '-c':
                shake_cmds[ix + 2] = '"%s"' % action
            else:
                shake_cmds.insert(ix + 1, '-c')
                shake_cmds.insert(ix + 2, '"%s"' % action)
            break

        cmd = config['shake_command'].replace('shake', config['shake_path'])
        cmd = cmd.replace('<EVID>', eventid)
        cmd = cmd.split() + shake_cmds
        p = subprocess.Popen(cmd)

    children[eventid] = {'popen': p, 'start_time': time.time()}

    return


def reap_children(children, config, logger):
    """
    Look through the list of child processes, reap the ones that have
    finished, and kill any that are taking too long.

    Args:
        children (dict): The dictionary of child processes and their
                         handles.
        config (dict): The configuration dictionary.
        logger (logger): The process logger.

    Returns:
        nothing: Nothing. Completed or killed child processes are removed
        from the list of children.
    """
    to_delete = []
    current_time = time.time()
    for eventid, info in children.items():
        returncode = info['popen'].poll()
        if returncode is not None:
            logger.info('Reaped child for event %s (return code %d)' %
                        (eventid, returncode))
            to_delete.append(eventid)
            continue
        #
        # Kill children who take too long
        #
        if info['start_time'] + config['max_process_time'] < current_time:
            logger.warning('Event %s taking too long, killing' % eventid)
            info['popen'].kill()
            info['popen'].wait()
            logger.warning('Reaped child for killed event %s' % eventid)
            to_delete.append(eventid)

    for eventid in to_delete:
        del children[eventid]

    return
