# System imports
import os
import os.path
import socket
import json
from collections import OrderedDict
import sys
from datetime import datetime, timezone
import time as time
import psutil
import copy
import shutil
import logging
from logging.handlers import TimedRotatingFileHandler
import subprocess
import shlex
import sqlite3

# Third-party imports
import daemon
import lockfile
from shapely.geometry import Point
from configobj import ConfigObj
from validate import Validator
from shapely.geometry import Polygon

# Local imports
from shakemap.utils.config import (get_config_paths,
                                   get_configspec,
                                   config_error)
from shakelib.rupture.origin import write_event_file
from shakemap.utils.amps import AmplitudeHandler
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

    network_delays = {}
    for key, value in config['network_delays'].items():
        network_delays[key] = str_to_seconds(value)
    config['network_delays'] = network_delays

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


class Queue(object):

    def __init__(self, pargs):

        current_time = int(time.time())
        self.MEMORY_UPDATE_TIME = current_time
        self.ASSOCIATE_UPDATE_TIME = current_time
        self.DB_MAINTENANCE_TIME = current_time

        self.children = {}
        self.attached = pargs.attached

        self.install_path, self.data_path = get_config_paths()

        self.config = get_config(self.install_path)
        #
        # Get shake.conf for the autorun modules
        #
        config_file = os.path.join(self.install_path, 'config', 'shake.conf')
        spec_file = get_configspec('shake')
        shake_config = ConfigObj(config_file, configspec=spec_file)
        results = shake_config.validate(Validator())
        if not isinstance(results, bool) or not results:
            config_error(shake_config, results)
        self.shake_cmds = shlex.split(shake_config['autorun_modules'])
        #
        # Turn this process into a daemon
        #
        self.logpath = os.path.join(self.install_path, 'logs')
        if not os.path.isdir(self.logpath):
            os.makedirs(self.logpath)
        pidfile = os.path.join(self.logpath, 'queue.pid')
        self.filelock = lockfile.FileLock(pidfile)
        if self.filelock.is_locked():
            if pargs.break_lock:
                self.filelock.break_lock()
            else:
                logger = self.getLogger()
                logger.error("pid lock file '%s' exists, can't start "
                             "sm_queue; exiting..." % (pidfile))
                sys.exit(-1)

    def queueMainLoop(self):

        context = daemon.DaemonContext(working_directory=self.data_path,
                                       pidfile=self.filelock)

        with self.getContext(context):
            self.logger = self.getLogger()
            #
            # Create the database for running and queued events.
            #
            self.eventQueue = EventQueue(self.install_path)
            #
            # Create the socket
            #
            qsocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            qsocket.bind(('', self.config['port']))
            # Set a timeout so that we can occasionally look for other
            # things to do
            qsocket.settimeout(30)
            qsocket.listen(5)
            #
            # Get a connection to the event database
            #
            self.ampHandler = AmplitudeHandler(self.install_path,
                                               self.data_path)

            self.logger.info('sm_queue initiated')

            #
            # At startup we want to see what was running when this process
            # shut down and try to restart them.
            #
            running = self.eventQueue.getRunningEvents()
            for eventid, command in running:
                self.logger.info("Startup: Running event %s" % (eventid))
                event = self.ampHandler.getEvent(eventid)
                # Update the XML because the DB may have newer information
                self.writeEventXml(event)
                p = subprocess.Popen(command)
                self.children[eventid] = {'popen': p,
                                          'start_time': time.time()}

            while True:
                #
                # Do routine stuff
                #
                self.doPeriodicTasks()
                #
                # Now wait for a connection
                #
                try:
                    (clientsocket, address) = qsocket.accept()
                except socket.timeout:
                    #
                    # Normal timeout; do routine tasks and then go
                    # back to waiting for a connection
                    #
                    continue
                #
                # Got a connection
                #
                hostname, _, _ = socket.gethostbyaddr(address[0])
    #            hostname = socket.getfqdn(hostname)
                self.logger.info('Got connection from %s at port %s' %
                                 (hostname, address[1]))

                if hostname not in self.config['servers']:
                    self.logger.warning('Connection from %s refused: not in '
                                        'valid servers list' % hostname)
                    clientsocket.close()
                    continue

                #
                # The accept() should guarantee that there's something
                # to read, but something could go wrong...
                #
                try:
                    clientsocket.settimeout(5)
                    data = clientsocket.recv(MAX_SIZE)
                except socket.timeout:
                    self.logger.warning('Did not get data from connection, '
                                        'continuing')
                    clientsocket.close()
                    continue
                else:
                    clientsocket.close()
                #
                # Decode the data and do something
                #
                try:
                    cmd = json.loads(data.decode('utf8'))
                except json.decoder.JSONDecodeError:
                    self.logger.warning("Couldn't decode data from %s: "
                                        "ignoring" % hostname)
                    continue

                if not isinstance(cmd, dict) or 'type' not in cmd or \
                   'data' not in cmd or 'id' not in cmd['data']:
                    self.logger.warning('Bad data from %s: ignoring' %
                                        hostname)
                    continue

                if cmd['type'] == 'origin':
                    self.logger.info('Received "origin" for event %s' %
                                     cmd['data']['id'])
                    if 'action' in cmd['data']:
                        action = cmd['data']['action']
                    else:
                        action = 'Origin received'
                    self.processOrigin(cmd['data'], action)
                elif cmd['type'] == 'cancel':
                    self.logger.info('Received "cancel" for event %s' %
                                     cmd['data']['id'])
                    self.processCancel(cmd['data'])
                else:
                    self.logger.info('Received "%s" for event %s' %
                                     cmd['type'], cmd['data']['id'])
                    self.processOther(cmd['data'], cmd['type'])

    def doPeriodicTasks(self):
        """ Check for finished children and start any needed timed repeats.

        Returns:
            nothing: Nothing.
        """
        #
        # Do routine stuff:
        # first check for repeats to queue
        #
        current_time = int(time.time())
        repeats = self.ampHandler.getRepeats()
        for eventid, otime, rep_list in repeats:
            while rep_list is not None and rep_list[0] < current_time:
                event = self.ampHandler.getEvent(eventid)
                if eventid in self.children:
                    #
                    # Event is already running; pop this repeat and move on
                    #
                    rep_list.pop(0)
                    if len(rep_list) == 0:
                        rep_list = None
                    event['repeats'] = rep_list
                    self.ampHandler.insertEvent(event, update=True)
                    continue
                rep_list.pop(0)
                if len(rep_list) == 0:
                    rep_list = None
                event['repeats'] = rep_list
                self.ampHandler.insertEvent(event, update=True)
                if event['lastrun'] == 0:
                    # This is a delayed first run
                    self.logger.info('Queueing event %s after network delay' %
                                     eventid)
                    self.dispatchEvent(event, 'Event added')
                else:
                    self.logger.info('Queueing repeat of event %s' % eventid)
                    self.dispatchEvent(event, 'Scheduled repeat')
                break

        #
        # Run the associator and dispatch events with new data
        #
        if self.config['associate_interval'] >= 0 and \
                self. ASSOCIATE_UPDATE_TIME + \
                self.config['associate_interval'] < current_time:
            self.ASSOCIATE_UPDATE_TIME = current_time
            self.associateAll()

        #
        # Reap any dead children, and then try to run queued events
        #
        _ = self.reapChildren()
        self.runQueuedEvents()

        #
        # Print memory usage once per hour to see how much we're leaking...
        #
        if self.MEMORY_UPDATE_TIME + 3600 < current_time:
            self.MEMORY_UPDATE_TIME = current_time
            process = psutil.Process(os.getpid())
            mem = getattr(process.memory_full_info(), 'uss', 0) / 1048576.0
            self.logger.info('Currently using %.1f MB' % mem)

        #
        # Do the occasional DB cleanup once per day; keep amps for 30
        # days and events for 1 year
        #
        if self.DB_MAINTENANCE_TIME + 86400 < current_time:
            self.DB_MAINTENANCE_TIME = current_time
            #
            # First do the assocication to make sure we don't drop any
            # amps that might associate
            #
            if self.config['associate_interval'] >= 0:
                self.ASSOCIATE_UPDATE_TIME = current_time
                self.associateAll()
            #
            # Now clean out the amps and events
            #
            self.ampHandler.cleanAmps(threshold=30)
            self.ampHandler.cleanEvents(threshold=365)

        return

    def getLogger(self):
        """Set up a logger for this process.

        Returns:
            logging.logger: An instance of a logger.
        """
        if not os.path.isdir(self.logpath):
            os.makedirs(self.logpath)
        logger = logging.getLogger('queue_logger')
        logger.setLevel(logging.INFO)
        if self.attached:
            handler = logging.StreamHandler()
        else:
            logfile = os.path.join(self.logpath, 'queue.log')
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

    def associateAll(self):
        """Do the associateAll method of the the AmplitudeHandler and
        then process all of the events with updated data.

        Returns:
            nothing: Nothing.
        """
        event_list = self.ampHandler.associateAll(pretty_print=True)
        for eventid in event_list:
            event = self.ampHandler.getEvent(eventid)
            self.processOrigin(event, 'Data association')

        return

    def writeEventXml(self, event):
        """ Create the event directory if it doesn't exist and write/re-write
        the event.xml file

        Args:
            event (dict): The event data structure.

        Returns:
            nothing: Nothing.
        """
        ttemp = event['time']
        try:
            dt = datetime.strptime(ttemp, constants.TIMEFMT)
        except ValueError:
            try:
                dt = datetime.strptime(ttemp, constants.ALT_TIMEFMT)
            except ValueError:
                self.logger.error("Can't parse input time %s" % ttemp)
                return
        event['time'] = dt

        event_dir = os.path.join(self.data_path, event['id'], 'current')
        if not os.path.isdir(event_dir):
            os.makedirs(event_dir)
        event_xml = os.path.join(event_dir, 'event.xml')

        self.logger.info('Writing event %s to event.xml' % (event['id']))
        val = write_event_file(event, event_xml)
        if val:
            self.logger.error('Error writing event.xml: %s' % val)

        event['time'] = ttemp

        return

    def moveEventDirectory(self, oldid, newid):
        """Change the name of an existing event directory to a new ID.
        """
        try:
            shutil.move(os.path.join(self.data_path, oldid),
                        os.path.join(self.data_path, newid))
        except shutil.Error as e:
            self.logger("Error trying to move data directory %s to %s: %s" %
                        (oldid, newid, str(e)))
        return

    def processOrigin(self, event, action):
        """ Determine if an event should be processed (or reprocessed) and
        dispatch it for processing.

        Args:
            event (dict): The event data structure.
            action (str): The "type" of the trigger that caused this function
                          to be called.

        Returns:
            nothing: Nothing.
        """
        current_time = int(time.time())
        force_run = False
        dispatch = True
        #
        # See if we already have this event, make a decision
        #
        existing = self.ampHandler.getEvent(event['id'])
        if existing is None and 'alt_eventids' in event:
            #
            # We haven't processed this ID, but the ID may have changed
            #
            for eid in event['alt_eventids'].split(','):
                if eid == event['id']:
                    continue
                alt_exists = self.ampHandler.getEvent(eid)
                if alt_exists is None:
                    continue
                #
                # We processed this event under a different ID
                # If the event is currently running with the old ID, kill it
                #
                if eid in self.children:
                    self.children[eid]['popen'].kill()
                    self.children[eid]['popen'].wait()
                    del self.children[eid]
                    self.eventQueue.deleteRunningEvent(eid)
                # Delete the old event from the database
                self.ampHandler.deleteEvent(eid)
                # Move the old event directory to the new ID
                self.moveEventDirectory(eid, event['id'])
                # Now treat the the new event ID as a new event.
                existing = None
                # But force it to run (because we want to update the event ID),
                # bypassing date and magnitude checks that come later
                force_run = True
                break

        if existing is None:
            #
            # This is a new event (or an event with a changed ID)
            #
            update = False
            # Do we want to run this event?
            if not force_run and self.magnitudeTooSmall(event):
                self.logger.info('Event %s (mag=%f) too small, skipping' %
                                 (event['id'], event['mag']))
                return
            if not force_run and self.eventTooOldOrInFuture(event):
                self.logger.info('Event %s too old or too far in the future, '
                                 'skipping' % event['id'])
                return
            #
            # Looks like we'll be running this event, get the repeats
            # (if any) and toss the ones that have already passed
            #
            replist = None
            try:
                dt = datetime.strptime(event['time'], constants.TIMEFMT)
            except ValueError:
                try:
                    dt = datetime.strptime(event['time'],
                                           constants.ALT_TIMEFMT)
                except ValueError:
                    self.logger.error("Can't parse input time %s" %
                                      event['time'])
                    return
            event_timestamp = int(dt.replace(tzinfo=timezone.utc).timestamp())
            for mag in sorted(self.config['repeats'].keys(), reverse=True):
                if event['mag'] > mag:
                    replist = [x + event_timestamp for x in
                               self.config['repeats'][mag]
                               if event_timestamp + x > current_time]
                    break
            #
            # The first time we run an event, we need to check its
            # network ID against those in the delay list. If present,
            # we add the required delay as the first repeat, but
            # don't dispatch the event. If the delay time has already
            # passed, just treat this as a normal event
            #
            if event['netid'] in self.config['network_delays']:
                delay = self.config['network_delays'][event['netid']]
                if event_timestamp + delay > current_time:
                    self.logger.info('Delaying processing event %s due to '
                                     'network delay configuration.' %
                                     (event['id']))
                    replist.insert(0, event_timestamp + delay)
                    dispatch = False

            event['repeats'] = replist if len(replist) > 0 else None
            event['lastrun'] = 0
        else:
            #
            # We've run this event before
            #
            update = True
            #
            # We want to update the event info in the database but
            # save the lastrun and repeats settings
            #
            event['lastrun'] = existing['lastrun']
            event['repeats'] = copy.copy(existing['repeats'])
        #
        # Insert or update the event info in the database, then
        # possibly queue the event to be run.
        #
        self.ampHandler.insertEvent(event, update=update)
        if dispatch is True:
            self.dispatchEvent(event, action)

        return

    def processOther(self, data, action):
        """A trigger has been issued for an event. Treat this as
        an origin update. If the event in question is not in our
        database, ignore the message.

        Args:
            data (dict): The event information dictionary.
            action (str): The "type" of the trigger that caused this function
                          to be called.

        Returns:
            nothing: Nothing.
        """
        eventid = data['id']
        existing = self.ampHandler.getEvent(eventid)
        if existing:
            self.processOrigin(existing, action)
        else:
            if 'alt_eventids' in data:
                for eid in data['alt_eventids'].split(','):
                    if eid == eventid:
                        continue
                    existing = self.ampHandler.getEvent(eid)
                    if existing:
                        self.processOrigin(existing, action)
                        return
            self.logger.info('Trigger of action "%s" is for unprocessed '
                             'event %s: ignoring' % (action, data['id']))
        return

    def processCancel(self, data):
        """We've received a cancellation of an event: run 'shake cancel'.

        Args:
            data (dict): The dictionary must have an event ID
                         under the 'id' key.

        Returns:
            nothing: Nothing.
        """
        eventid = data['id']
        existing = self.ampHandler.getEvent(eventid)
        if existing:
            self.dispatchEvent(data, 'cancel')
            return

        if 'alt_eventids' in data:
            for eid in data['alt_eventids'].split(','):
                if eid == eventid:
                    continue
                existing = self.ampHandler.getEvent(eid)
                if existing:
                    self.dispatchEvent(existing, 'cancel')
                    return

        self.logger.info('cancel is for unprocessed event %s: ignoring' %
                         eventid)
        return

    def magnitudeTooSmall(self, event):
        """ Return False if the magnitude is greater than the threshold
        magnitude of the first metro box within which it falls, or
        the global minmag if it does not fall within a box; return
        true otherwise.

        Args:
            event (dict): The event dictionary; must contain at least "mag",
                "lon", and "lat" keys.

        Returns:
            bool: True if the event is too small to process; False otherwise.
        """
        mag = event['mag']
        lon = event['lon']
        lat = event['lat']
        pt = Point((lon, lat))
        for boxname in sorted(self.config['boxes']):
            boxdict = self.config['boxes'][boxname]
            if pt.within(boxdict['poly']):
                if mag >= boxdict['mag']:
                    return False
                else:
                    return True
        #
        # Not in any boxes
        #
        if mag >= self.config['minmag']:
            return False

        return True

    def eventTooOldOrInFuture(self, event):
        """ Return True if the event is too old or too far in the future to
        process; return False otherwise.

        Args:
            event (dict): The event data structure.

        Returns:
            bool: True if the event is older than old_event_age or is
            more than future_event_age in the future; returns False otherwise.
        """

        current_time = time.time()
        try:
            event_time = datetime.strptime(event['time'], constants.TIMEFMT).\
                replace(tzinfo=timezone.utc).timestamp()
        except ValueError:
            event_time = datetime.strptime(event['time'],
                                           constants.ALT_TIMEFMT).\
                replace(tzinfo=timezone.utc).timestamp()
        if self.config['old_event_age'] >= 0 and \
           event_time + self.config['old_event_age'] < current_time:
            return True
        if self.config['future_event_age'] >= 0 and \
           event_time - self.config['future_event_age'] > current_time:
            return True

        return False

    def dispatchEvent(self, event, action):
        """ Queue a run for the specified event.

        Args:
            event (dict): The data structure of the event to process.
            action (str): 'cancel', 'test', or some other string. 'cancel'
                          starts the cancel process, 'test' queues the process
                          'echo eventid'. Any other string queues the
                          shake process to be run at the next opportunity.
                          See the configuration file 'queue.conf' for the
                          exact commands that will be run.

        Returns:
            nothing: Nothing.
        """
        eventid = event['id']
        if action == 'cancel':
            #
            # Cancellations aren't queued, they're run immediately
            #
            self.logger.info('Canceling event %s' % eventid)
            if eventid in self.children:
                self.logger.info('Event %s is running; killing...' % eventid)
                self.children[eventid]['popen'].kill()
                self.children[eventid]['popen'].wait()
                del self.children[eventid]
                self.eventQueue.deleteRunningEvent(eventid)
            cmd = self.config['cancel_command'].replace(
                'shake', self.config['shake_path'])
            cmd = cmd.replace('<EVID>', eventid)
            cmd = cmd.split()
            p = subprocess.Popen(cmd)
            self.children[eventid] = {'popen': p, 'start_time': time.time()}
            self.eventQueue.insertRunningEvent(eventid, cmd)
            return

        self.logger.info('Queueing event %s due to action "%s"' %
                         (eventid, action))
        #
        # Add the action as the assemble/augment comment, or replace the
        # comment if it is already there.
        #
        for ix, shcmd in enumerate(self.shake_cmds):
            if shcmd not in ['assemble', 'augment']:
                continue
            if len(self.shake_cmds) == ix + 1:  # This shouldn't happen
                self.shake_cmds.append('-c')
                self.shake_cmds.append('"%s"' % action)
            elif self.shake_cmds[ix + 1] == '-c':
                self.shake_cmds[ix + 2] = '"%s"' % action
            else:
                self.shake_cmds.insert(ix + 1, '-c')
                self.shake_cmds.insert(ix + 2, '"%s"' % action)
            break

        if action == 'test':
            cmd = self.config['shake_command'].replace('shake', 'echo')
        else:
            cmd = self.config['shake_command'].replace(
                'shake', self.config['shake_path'])
        cmd = cmd.replace('<EVID>', eventid)
        cmd = cmd.split() + self.shake_cmds

        self.eventQueue.queueEvent(eventid, cmd, event['mag'])
        return

    def runQueuedEvents(self):
        """If there is space, run events from the queue
        """
        if len(self.children) >= self.config['max_subprocesses']:
            self.logger.info('Processing queue is full; waiting for open '
                             'slots.')
            return
        current_time = int(time.time())
        mtw = self.config['max_trigger_wait']
        queued = self.eventQueue.getQueuedEvents()
        for eventid, command in queued:
            event = self.ampHandler.getEvent(eventid)
            if eventid in self.children:
                #
                # Event is currently running, don't run it but make sure
                # there's a repeat pretty soon
                #
                if event['repeats']:
                    if event['repeats'][0] > current_time + mtw:
                        event['repeats'].insert(0, current_time + mtw)
                else:
                    event['repeats'] = [current_time + mtw]
                self.ampHandler.insertEvent(event, update=True)
                self.logger.info('Event %s is currently running, shelving '
                                 'this update' % event['id'])
                self.eventQueue.dequeueEvent(eventid)
                continue
            if event['repeats']:
                delta_t = current_time - event['repeats'][0]
                if delta_t > -mtw:
                    # We're due for a rerun anyway, so just leave the
                    # event queued
                    self.logger.info('Event %s will repeat soon, shelving '
                                     'this update' % eventid)
                    self.eventQueue.dequeueEvent(eventid)
                    continue
            if current_time - event['lastrun'] < mtw:
                #
                # We ran this event very recently, but don't have a repeat
                # scheduled in the near future, so let's skip this one
                # but make sure something happens relatively soon
                #
                if event['repeats']:
                    event['repeats'].insert(0, current_time + mtw)
                else:
                    event['repeats'] = [current_time + mtw]
                self.ampHandler.insertEvent(event, update=True)
                self.logger.info('Event %s ran recently, shelving this '
                                 'update' % event['id'])
                self.eventQueue.dequeueEvent(eventid)
                continue

            self.logger.info("Running event %s" % (eventid))
            # Update the XML because the DB may have newer information
            self.writeEventXml(event)
            p = subprocess.Popen(command)
            self.children[eventid] = {'popen': p, 'start_time': time.time()}
            self.eventQueue.dequeueEvent(eventid)
            self.eventQueue.insertRunningEvent(eventid, command)
            if len(self.children) >= self.config['max_subprocesses']:
                self.logger.info('Processing queue is full; waiting for open '
                                 'slots.')
                break

        return

    def reapChildren(self):
        """
        Look through the list of child processes, reap the ones that have
        finished, and kill any that are taking too long.

        Returns:
            nothing: Nothing. Completed or killed child processes are removed
            from the list of children.
        """
        to_delete = []
        current_time = time.time()
        for eventid, info in self.children.items():
            returncode = info['popen'].poll()
            if returncode is not None:
                self.logger.info('Reaped child for event %s (return code %d)' %
                                 (eventid, returncode))
                event = self.ampHandler.getEvent(eventid)
                if event:
                    event['lastrun'] = current_time
                    self.ampHandler.insertEvent(event, update=True)
                to_delete.append(eventid)
                self.eventQueue.deleteRunningEvent(eventid)
                continue
            #
            # Kill children who take too long
            #
            if info['start_time'] + self.config['max_process_time'] < \
                    current_time:
                self.logger.warning('Event %s taking too long, killing' %
                                    eventid)
                info['popen'].kill()
                info['popen'].wait()
                self.logger.warning('Reaped child for killed event %s' %
                                    eventid)
                to_delete.append(eventid)
                self.eventQueue.deleteRunningEvent(eventid)

        for eventid in to_delete:
            del self.children[eventid]

        return to_delete

    def getContext(self, context):
        """Returns a context based on the value of the 'attached' argument.
        If attached is True, then the function returns an instance of the
        Dummycontext; if it is False the function returns the 'context'
        argument.

        Args:
            context (Context manager): A valid context manager.

        Returns:
            Context manager: If attached is True, the function returns an
                             instance of the Dummycontext; if False, returns
                             the 'context' argument.
        """
        if self.attached:
            return Dummycontext()
        else:
            return context


class Dummycontext(object):
    """This is a dummy context that can be used as a context manager
    with 'with'. It doesn't do anything.
    """
    def __enter__(self): return self

    def __exit__(*x): pass


class EventQueue(object):
    """Class to maintain some persistence in the event of a program
    crash or shutdown. The db file can be removed if the operator wants
    a fresh start.
    """
    def __init__(self, ipath):
        queued_events = OrderedDict([('id', 'INTEGER_PRIMARY_KEY'),
                                     ('eventid', 'TEXT UNIQUE'),
                                     ('command', 'TEXT'),
                                     ('mag', 'REAL')])
        running_events = OrderedDict([('id', 'INTEGER_PRIMARY_KEY'),
                                      ('eventid', 'TEXT UNIQUE'),
                                      ('command', 'TEXT')])
        tables = {'queued': queued_events,
                  'running': running_events}
        self.db_file = os.path.join(ipath, 'data', 'event_queue.db')
        db_exists = os.path.isfile(self.db_file)
        self._connection = sqlite3.connect(self.db_file, timeout=15)
        if self._connection is None:
            raise RuntimeError('Could not connect to %s' % self.db_file)
        self._connection.isolation_level = 'EXCLUSIVE'
        self._cursor = self._connection.cursor()
        self._cursor.execute('PRAGMA foreign_keys = ON')
        self._cursor.execute('PRAGMA journal_mode = WAL')
        if not db_exists:
            for table, tdict in tables.items():
                createcmd = 'CREATE TABLE %s (' % table
                nuggets = []
                for column, ctype in tdict.items():
                    nuggets.append('%s %s' % (column, ctype))
                createcmd += ','.join(nuggets) + ')'
                self._cursor.execute(createcmd)
            self._cursor.execute('CREATE INDEX queue_index ON '
                                 'queued(eventid)')
            self._cursor.execute('CREATE INDEX mag_index ON '
                                 'queued(mag)')
            self._cursor.execute('CREATE INDEX event_index ON '
                                 'running(eventid)')
            self._cursor.execute('PRAGMA journal_mode = WAL')

    def __del__(self):
        """Destructor.

        """
        if hasattr(self, '_connection') and self._connection is not None:
            self._disconnect()

    def _disconnect(self):
        self.commit()
        self._cursor.close()
        self._connection.close()
        self._connection = None
        self._cursor = None

    def commit(self):
        """Commit any operations to the database.
        """
        self._connection.commit()

    def getQueuedEvents(self):
        query = 'SELECT eventid, command, mag FROM queued ORDER BY mag DESC'
        self._cursor.execute(query)
        erows = self._cursor.fetchall()
        return [(x[0], json.loads(x[1])) for x in erows]

    def queueEvent(self, eventid, command, mag):
        query = 'REPLACE INTO queued (eventid, command, mag) VALUES (?, ?, ?)'
        self._cursor.execute(query, (eventid, json.dumps(command), mag))
        self.commit()

    def dequeueEvent(self, eventid):
        query = 'DELETE FROM queued WHERE eventid = ?'
        self._cursor.execute(query, (eventid,))
        self.commit()

    def getRunningEvents(self):
        query = 'SELECT eventid, command FROM running'
        self._cursor.execute(query)
        erows = self._cursor.fetchall()
        return [(x[0], json.loads(x[1])) for x in erows]

    def insertRunningEvent(self, eventid, command):
        query = 'INSERT INTO running (eventid, command) VALUES (?, ?)'
        self._cursor.execute(query, (eventid, json.dumps(command)))
        self.commit()

    def deleteRunningEvent(self, eventid):
        query = 'DELETE FROM running WHERE eventid = ?'
        self._cursor.execute(query, (eventid,))
        self.commit()
