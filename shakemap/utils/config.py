# stdlib imports
import os
import os.path
import pkg_resources
import logging
import logging.config

# third party libraries
import numpy as np
from configobj import ConfigObj, flatten_errors
from validate import Validator, ValidateError

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


def get_data_path():
    """
    Return the path to the shakemap package data directory holding
    configspec.conf, the template config files, and the example config
    files.

    Returns:
        str: The full path to the data directory.

    """
    return pkg_resources.resource_filename('shakemap', 'data')


def get_configspec(config=None):
    """
    Returns the full path to a ShakeMap config spec file.

    Args:
      config (str): Name of config spec to find, or None.

    Returns:
        str: The path to a config spec, or

    """
    if config is None:
        return os.path.join(get_data_path(), 'configspec.conf')
    fname = os.path.join(get_data_path(), '%sspec.conf' % config)
    if not os.path.isfile(fname):
        return FileNotFoundError('No file "%s" exists.' % fname)
    return fname


def get_config_paths():
    """
    Returns two paths based on the currently selected profile in the
    user's ~/.shakemap/profile.conf: 1) the path to the ShakeMap
    installation directory; 2) the path to the data directory.

    If this function is called within a pytest process, it will
    return the paths to the repository's test install and data
    directories.

    Returns:
        tuple: The paths to the ShakeMap install directory
        and the data directory.
    """
    if 'CALLED_FROM_PYTEST' in os.environ:
        base_path = os.path.join(get_data_path(), '..', '..', 'tests', 'data')
        install = os.path.join(base_path, 'install')
        data = os.path.join(base_path, 'eventdata')
    else:
        config_file = os.path.join(os.path.expanduser('~'), '.shakemap',
                                   'profiles.conf')
        config = ConfigObj(config_file)
        profile_name = config['profile']
        profile = config['profiles'][profile_name]
        install = profile['install_path']
        data = profile['data_path']
    return (install, data)


def get_custom_validator():
    """
    Returns a validator suitable for validating the ShakeMap config
    files.

    Returns:
        :class:`Validator`: A Validator object.

    """
    fdict = {
        'file_type': file_type,
        'directory_type': directory_type,
        'annotatedfloat_type': annotatedfloat_type,
        'gmpe_list': gmpe_list,
        'weight_list': weight_list,
        'extent_list': extent_list,
        'status_string': status_string,
    }
    validator = Validator(fdict)
    return validator


def config_error(config, results):
    """
    Parse the results of a ConfigObj validation and print the errors.
    Throws a RuntimeError exception  upon completion if any errors or
    missing sections are encountered.

    Args:
        config (ConfigObj): The ConfigObj instance representing the
            parsed config.
        results (dict): The dictionary returned by the validation of
            the 'config' arguments.

    Returns:
        Nothing: Nothing

    Raises:
        RuntimeError: Should always raise this exception.
    """
    errs = 0
    for (section_list, key, _) in flatten_errors(config, results):
        if key is not None:
            print('The "%s" key in the section "%s" failed validation'
                  % (key, ', '.join(section_list)))
            errs += 1
        else:
            print('The following section was missing:%s '
                  % ', '.join(section_list))
            errs += 1
    if errs:
        raise RuntimeError('There %s %d %s in configuration.'
                           % ('was' if errs == 1 else 'were', errs,
                              'error' if errs == 1 else 'errors'))


def check_config(config, logger):
    """
    Checks that the gmpe, gmice, ipe, ccf, and component parameters
    in config are defined in their respective sections. Raises a
    ValidateError exception if an error is encountered.

    Args:
        config (ConfigObj): A ConfigObj instance.
        logger (logger): The logger to which to write complaints.

    Returns:
        Nothing: Nothing.

    """
    if config['modeling']['gmpe'] not in config['gmpe_sets']:
        logger.error('Configuration error: gmpe %s not in gmpe_sets' %
                     (config['modeling']['gmpe']))
        raise ValidateError()
    if config['modeling']['gmice'] not in config['gmice_modules']:
        logger.error('Configuration error: gmice %s not in gmice_modules' %
                     (config['modeling']['gmice']))
        raise ValidateError()
    if config['modeling']['ipe'] not in config['ipe_modules']:
        logger.error('Configuration error: ipe %s not in ipe_modules' %
                     (config['modeling']['ipe']))
        raise ValidateError()
    if config['modeling']['ccf'] not in config['ccf_modules']:
        logger.error('Configuration error: ccf %s not in ccf_modules' %
                     (config['modeling']['ccf']))
        raise ValidateError()
    if config['interp']['component'] not in config['component_modules']:
        logger.error('Configuration error: component %s not in '
                     'component_modules' %
                     (config['interp']['component']))
        raise ValidateError()


def annotatedfloat_type(value):
    """
    Checks to see if value is a float, or a float with a 'c', 'm', or 'd'
    appended. Then converts the value to decimal degrees where an unadorned
    float or a float plus 'd' is interpreted as decimal degrees, 'm' is
    interpreted as arc-minutes, and 'c' is interpreted as arc-seconds.
    Raises a ValidateError exception on failure.

    Args:
        value (str): A string representing a float or a float appended
            with 'd', 'm', or 'c' (for degrees, minutes, seconds).

    Returns:
        float: The input value converted to decimal degrees.

    """
    try:
        out = float(value)
    except ValueError:
        try:
            if value.endswith('c'):
                out = float(value.replace('c', '')) / 3600.0
            elif value.endswith('m'):
                out = float(value.replace('m', '')) / 60.0
            elif value.endswith('d'):
                out = float(value.replace('d', ''))
            else:
                raise ValidateError(value)
        except:
            raise ValidateError(value)
    return out


def weight_list(value, min):
    """
    Checks to see if value is a list of floats at least min elements long,
    and whose values add up to 1.0.  Raises a ValidateError exception on
    failure.

    Args:
        value (str): A string representing a list of floats.

    Returns:
        list: The input string converted to a list of floats.

    """

    if isinstance(value, str) and value == 'None':
        if int(min) == 0:
            return []
        else:
            print("list must contain at least %d entr%s" %
                  (min, "ies" if int(min) != 1 else "y"))
            raise ValidateError()
    if isinstance(value, str):
        if value.startswith('[') and value.endswith(']'):
            value = value.replace('[', '')
            value = value.replace(']', '')
        if not value:
            if int(min) == 0:
                value = []
            else:
                print("list must contain at least %d entr%s" %
                      (min, "ies" if int(min) != 1 else "y"))
                raise ValidateError()
        else:
            value = [value]
    if len(value) < int(min):
        print("list must contain at least %d entr%s" %
              (min, "ies" if int(min) != 1 else "y"))
        raise ValidateError()
    try:
        out = [float(a) for a in value]
    except:
        print("%s is not a list of floats" % value)
        raise ValidateError()
    np_out = np.array(out)
    if np.any(np_out < 0):
        print("all list values must be >= 0: %s" % value)
        raise ValidateError()
    if len(out) > 0 and np.abs(np.sum(np_out) - 1.0) > 0.01:
        print("weights must sum to 1.0: %s" % value)
        raise ValidateError()

    return out


def gmpe_list(value, min):
    """
    Checks to see if value is a list of strings at least min elements long.
    The entries are not checked for their validity as GMPEs. Raises a
    ValidateError exception on failure.

    Args:
        value (str): A string representing a list of GMPEs.

    Returns:
        list: The input string converted to a list of GMPEs.

    """

    if value == 'None' or value == '[]':
        value = []
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, list) or len(value) < int(min):
        print("'%s' is not a list of at least %s gmpes" % (value, min))
        raise ValidateError()
    for gmpe in value:
        if not isinstance(gmpe, str):
            print("'%s' is not a list of strings" % (value))
            raise ValidateError()

    return value


def extent_list(value):
    """
    Checks to see if value is an empty list or a list of four floats,
    whose values are valid coordinates in (longitude, longitude,
    latitude, latitude) order. Returns a list upon success; raises a
    ValidateError exception on failure.

    Args:
        value (str): A string representing a list of geographic
            coordinates.

    Returns:
        list: The input string converted to a list of floats.

    """

    if isinstance(value, str) and (value == 'None' or value == '[]'):
        return []
    if isinstance(value, list) and not value:
        return []
    if not isinstance(value, list):
        print("'%s' is not a list of 4 coordinates" % value)
        raise ValidateError()
    if len(value) != 4:
        print("extent list must contain 4 entries")
        raise ValidateError()
    try:
        out = [float(a) for a in value]
    except:
        print("%s is not a list of 4 floats" % value)
        raise ValidateError()
    if out[0] < -360.0 or out[0] > 360.0 or \
       out[2] < -360.0 or out[2] > 360.0 or \
       out[1] < -90.0 or out[1] > 90.0 or \
       out[3] < -90.0 or out[3] > 90.0:
        print("Invalid extent: ", value,
              " : -360 <= longitude <= 360, -90 <= latitude <= 90")
        raise ValidateError()

    return out


def file_type(value):
    """
    Checks to see if value is a valid file or an empty string.
    Raises a ValidateError exception on failure.

    Args:
        value (str): A string representing the path to a file.

    Returns:
        str: The input string.

    """
    if not value or value == 'None':
        return ''
    if not os.path.isfile(value):
        print("file '%s' is not a valid file" % value)
        raise ValidateError(value)
    return value


def directory_type(value):
    """
    Checks to see if value is a valid directory or an empty string.
    Raises a ValidateError exception on failure.

    Args:
        value (str): A string representing the path to a directory.

    Returns:
        str: The input string.

    """
    if not value or value == 'None':
        return ''
    if not os.path.isdir(value):
        raise ValidateError(value)
    return value


def status_string(value, min):
    """
    Checks to see if value is one of the ShakeMap status string of
    'automatic', 'released', or 'reviewed.  Raises a ValidateError
    exception on failure.

    Args:
        value (str): A status string.

    Returns:
        str: The input string. 'automatic' is returned if value is empty.

    """
    if not value:
        return 'automatic'
    if value not in ('automatic', 'released', 'reviewed'):
        raise ValidateError(value)
    return value


def cfg_float_list(value):
    """
    Converts (if possible) the input list (or string) to a list
    of floats. Raises ValidateError if the input can't be
    converted to a list of floats.

    Args:
        value (str or list): A string or list of strings to be
            converted to a list of floats.

    Returns:
        list: The input converted to a list of floats.

    Raises:
        ValidateError
    """
    if not value or value == 'None':
        print("'%s' is not a list of at least 1 float" % (value))
        raise ValidateError()
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, list) or len(value) < 1:
        print("'%s' is not a list of at least 1 float" % (value))
        raise ValidateError()
    fvalue = []
    for val in value:
        try:
            fval = float(val)
        except ValueError:
            print("'%s' is not a list of floats" % (value))
            raise ValidateError()
        fvalue.append(fval)
    return fvalue


def cfg_float(value):
    """
    Converts (if possible) the input string to a float. Raises
    ValidateError if the input can't be converted to a float.

    Args:
        value (str): A string to be converted to a float.

    Returns:
        float: The input converted to a float.

    Raises:
        ValidateError
    """
    if not isinstance(value, (str, float)) or not value or value == 'None':
        print("'%s' is not a float" % (value))
        raise ValidateError()
    try:
        fval = float(value)
    except ValueError:
        print("'%s' is not a float" % (value))
        raise ValidateError()
    return fval


def get_shake_config():
    """
    Return a dictionary containing input required to run the shake program,
    including logging.

    See https://docs.python.org/3.5/library/logging.config.html.

    Returns:
       dict: Dictionary containing shake config information.
    """
    install_path, data_path = get_config_paths()
    conf_file = os.path.join(install_path, 'config', 'shake.conf')
    spec_file = get_configspec(config='shake')
    shake_conf = ConfigObj(conf_file,
                           configspec=spec_file,
                           interpolation='template')

    val = Validator()
    results = shake_conf.validate(val)
    if not isinstance(results, bool) or not results:
        config_error(global_config, results)

    return shake_conf


def get_logger(eventid, log_option=None):
    """Return the logger instance for ShakeMap.  Only use once!

    Args:
        eventid (str): Event ID.
        log_option (str): One of 'log','quiet', 'debug', or None.

    Returns:
        Logger: logging Logger instance.
    """
    install_path, data_path = get_config_paths()
    config = get_logging_config()
    if log_option == 'debug' or log_option == 'quiet' or log_option is None:
        format = config['formatters']['standard']['format']
        datefmt = config['formatters']['standard']['datefmt']
        # create a console handler, with verbosity setting chosen by user
        if log_option == 'debug':
            level = logging.DEBUG
        elif log_option == 'quiet':
            level = logging.ERROR
        elif log_option is None:  # default interactive
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
        global_log_dir = os.path.join(install_path, 'logs')
        if not os.path.isdir(global_log_dir):
            os.makedirs(global_log_dir, exist_ok=True)
        global_log_file = os.path.join(global_log_dir, 'shake.log')
        config['handlers']['global_file']['filename'] = global_log_file
        logging.config.dictConfig(config)
        log_cfg = list(config['loggers'])
    # get the root logger, otherwise we can't log in sub-libraries
    logger = logging.getLogger()

    return logger


def get_logging_config():
    """Extract logging configuration from shake config.

    See this URL for example of config.
    https://gist.github.com/st4lk/6287746

    See https://docs.python.org/3.5/library/logging.config.html

    Returns:
        dict: Dictionary suitable for use with logging.config.dictConfig().
    """

    shake_conf = get_shake_config()
    log_config = shake_conf['shake']
    _clean_log_dict(log_config)

    # Here follows a bit of trickery...
    # To have a logger point to the root logger using the dictConfig() method,
    # you need to have the logger have a name equal to the empty string ''.
    # Our logging dictionary is originally specified using ConfigObj, which
    # does not allow for empty section headers.  So, we need to get all of the
    # information from the logger we specify, copy it into a logger dictionary with
    # an empty key, and then delete the original logger from the config dictionary.
    # Whew.
    log_name = log_config['loggers'].keys()[0]
    log_config['loggers'][''] = log_config['loggers'][log_name]
    del log_config['loggers'][log_name]
    return log_config


def _clean_log_dict(config):
    """Clean up dictionary returned by ConfigObj into form suitable for logging.

    Basically, ConfigObj.validate wants all sections that are Handlers (for example)
    to have the same fields, so it fills them in with default values.  However,
    if you try to give a StreamHandler a filename parameter, it generates an error,
    hence the code below.

    Returns:
        dict: Dictionary suitable for use with logging.config.dictConfig().

    """
    for handlerkey, handler in config['handlers'].items():
        myclass = handler['class']
        req_fields = REQ_FIELDS[myclass]
        for key, value in handler.items():
            if key not in req_fields:
                del handler[key]
