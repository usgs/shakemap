# stdlib imports
import os
import os.path
import logging
import logging.config

# third party libraries
from configobj import ConfigObj
from validate import Validator

# local libraries
from shakemap.utils.config import (get_config_paths,
                                   get_configspec,
                                   config_error)

REQ_FIELDS = {
    'logging.handlers.TimedRotatingFileHandler': [
        'level',
        'formatter',
        'class',
        'when',
        'filename'
    ],
    'logging.FileHandler': [
        'level',
        'formatter',
        'class',
        'filename'],
    'logging.handlers.SMTPHandler': [
        'level',
        'formatter',
        'mailhost',
        'fromaddr',
        'toaddrs',
        'subject',
        'class']
}

# mapping of string versions of logging module logging levels
# to the corresponding constants.
LOG_LEVELS = {
    'DEBUG': logging.DEBUG,
    'INFO': logging.INFO,
    'WARNING': logging.WARNING,
    'ERROR': logging.ERROR,
    'CRITICAL': logging.CRITICAL
}


def get_logger(eventid, log_option=None, log_file=None):
    """Return the logger instance for ShakeMap.  Only use once!

    Args:
        eventid (str): Event ID.
        log_option (str): One of 'log','quiet', 'debug', or None.

    Returns:
        Logger: logging Logger instance.
    """
    install_path, data_path = get_config_paths()
    config = get_logging_config()
    if log_file is None:
        format = config['formatters']['standard']['format']
        datefmt = config['formatters']['standard']['datefmt']
        # create a console handler, with verbosity setting chosen by user
        if log_option == 'debug':
            level = logging.DEBUG
        elif log_option == 'quiet':
            level = logging.ERROR
        else:  # default interactive
            level = logging.INFO

        logdict = {
            'version': 1,
            'formatters': {
                'standard': {
                    'format': format,
                    'datefmt': datefmt
                }
            },
            'handlers': {
                'stream': {
                    'level': level,
                    'formatter': 'standard',
                    'class': 'logging.StreamHandler'
                }
            },
            'loggers': {
                '': {
                    'handlers': ['stream'],
                    'level': level,
                    'propagate': True
                }
            }
        }

        logging.config.dictConfig(logdict)
    else:
        event_log_dir = os.path.join(data_path, eventid)
        if not os.path.isdir(event_log_dir):
            raise NotADirectoryError("Can't open log file: event %s "
                                     "not found" % eventid)
        event_log_file = os.path.join(event_log_dir, 'shake.log')
        config['handlers']['event_file']['filename'] = event_log_file
        if log_option == 'debug':
            config['handlers']['event_file']['level'] = logging.DEBUG
        elif log_option == 'quiet':
            config['handlers']['event_file']['level'] = logging.ERROR
        global_log_dir = os.path.join(install_path, 'logs')
        if not os.path.isdir(global_log_dir):
            os.makedirs(global_log_dir, exist_ok=True)
        global_log_file = os.path.join(global_log_dir, 'shake.log')
        config['handlers']['global_file']['filename'] = global_log_file
        if log_option == 'debug':
            config['handlers']['global_file']['level'] = logging.DEBUG
        elif log_option == 'quiet':
            config['handlers']['global_file']['level'] = logging.ERROR
        if log_option == 'debug':
            config['loggers']['']['level'] = logging.DEBUG
        elif log_option == 'quiet':
            config['loggers']['']['level'] = logging.ERROR

        logging.config.dictConfig(config)

    # Have the logger capture anything from the 'warnings' package,
    # which many libraries use.
    logging.captureWarnings(True)

    # Get the root logger, otherwise we can't log in sub-libraries
    logger = logging.getLogger()

    return logger


def get_logging_config():
    """Extract logging configuration from logging.conf.

    See this URL for example of config.
    https://gist.github.com/st4lk/6287746

    See https://docs.python.org/3.5/library/logging.config.html

    Returns:
        dict: Dictionary suitable for use with logging.config.dictConfig().
    """

    install_path, data_path = get_config_paths()
    conf_file = os.path.join(install_path, 'config', 'logging.conf')
    spec_file = get_configspec(config='logging')
    log_config = ConfigObj(conf_file,
                           configspec=spec_file,
                           interpolation='template')

    val = Validator()
    results = log_config.validate(val)
    if not isinstance(results, bool) or not results:
        config_error(log_config, results)

    _clean_log_dict(log_config)

    # Here follows a bit of trickery...
    # To have a logger point to the root logger using the dictConfig() method,
    # you need to have the logger have a name equal to the empty string ''.
    # Our logging dictionary is originally specified using ConfigObj, which
    # does not allow for empty section headers.  So, we need to get all of the
    # information from the logger we specify, copy it into a logger dictionary
    # with an empty key, and then delete the original logger from the config
    # dictionary. Whew.
    log_name = log_config['loggers'].keys()[0]
    log_config['loggers'][''] = log_config['loggers'][log_name]
    del log_config['loggers'][log_name]
    return log_config


def _clean_log_dict(config):
    """Clean up dictionary returned by ConfigObj into form suitable for
    logging.

    Basically, ConfigObj.validate wants all sections that are Handlers (for
    example) to have the same fields, so it fills them in with default values.
    However, if you try to give a StreamHandler a filename parameter, it
    generates an error, hence the code below.

    Returns:
        dict: Dictionary suitable for use with logging.config.dictConfig().

    """
    for handlerkey, handler in config['handlers'].items():
        myclass = handler['class']
        req_fields = REQ_FIELDS[myclass]
        for key, value in handler.items():
            if key not in req_fields:
                del handler[key]
