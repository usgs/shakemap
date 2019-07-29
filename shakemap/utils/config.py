# stdlib imports
import os
import os.path
import sys
import pkg_resources
import logging
import glob

# third party libraries
import numpy as np
from configobj import (ConfigObj,
                       flatten_errors,
                       get_extra_values)
from validate import (Validator,
                      ValidateError)


def get_data_path():
    """
    Return the path to the shakemap package data directory holding
    modelspec.conf, the template config files, and the example config
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
        return os.path.join(get_data_path(), 'modelspec.conf')
    fname = os.path.join(get_data_path(), '%sspec.conf' % config)
    if not os.path.isfile(fname):
        raise FileNotFoundError('No file "%s" exists.' % fname)
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

    Raises:
        FileNotFoundError -- if the profile file can't be found.
        ValueError -- If the correct profile can't be foun in profiles.conf.
    """
    if 'CALLED_FROM_PYTEST' in os.environ:
        base_path = os.path.join(get_data_path(), '..', '..', 'tests', 'data')
        install = os.path.join(base_path, 'install')
        data = os.path.join(base_path, 'eventdata')
    else:
        config_file = os.path.join(
            os.path.expanduser('~'),
            '.shakemap',
            'profiles.conf')

        if not os.path.isfile(config_file):
            raise FileNotFoundError("Can't find a profile file: "
                                    "have you run sm_profile?")
        config = ConfigObj(config_file, encoding='utf-8-sig')

        config = check_profile_config(config)

        profile_name = config['profile']
        if not profile_name or profile_name == 'None':
            raise ValueError("No profile set in the profiles.conf file")
        if not config['profiles'] or profile_name not in config['profiles']:
            raise ValueError("Profile %s not found in list of profiles")
        profile = config['profiles'][profile_name]
        install = profile['install_path']
        data = profile['data_path']
    return (install, data)


def get_model_config(install_path, datadir, logger):
    #
    # Look for global configs in install_path/config
    #
    spec_file = get_configspec()
    validator = get_custom_validator()
    logger.debug('Looking for configuration files...')
    modules = ConfigObj(
        os.path.join(install_path, 'config', 'modules.conf'),
        configspec=spec_file)
    gmpe_sets = ConfigObj(
        os.path.join(install_path, 'config', 'gmpe_sets.conf'),
        configspec=spec_file)
    global_config = ConfigObj(
        os.path.join(install_path, 'config', 'model.conf'),
        configspec=spec_file)

    #
    # this is the event specific model.conf (may not be present)
    # prefer model.conf to model_select.conf
    #
    event_config_file = os.path.join(datadir, 'model.conf')
    event_config_zc_file = os.path.join(datadir, 'model_select.conf')
    if os.path.isfile(event_config_file):
        event_config = ConfigObj(event_config_file,
                                 configspec=spec_file)
    elif os.path.isfile(event_config_zc_file):
        event_config = ConfigObj(event_config_zc_file,
                                 configspec=spec_file)
    else:
        event_config = ConfigObj()

    #
    # start merging event_config
    #
    global_config.merge(event_config)
    global_config.merge(modules)
    global_config.merge(gmpe_sets)

    results = global_config.validate(validator)
    if not isinstance(results, bool) or not results:
        config_error(global_config, results)

    check_config(global_config, logger)

    return global_config


def path_macro_sub(s, ip='', dp='', gp='', ei=''):
    """
    Replace macros with current paths:

    - <INSTALL_DIR> is replaced with the contents of ip
    - <DATA_DIR> is replaced with the contents of dp
    - <GLOBAL_DATA> is replaced with the contents of gp
    - <EVENT_ID> is replaced with the contents of ei

    e.g., path_macro_sub("<INSTALL_DIR>/<DATA_DIR>", "hello", "world")
    would return "hello/world". It is not an error if the original string
    does not contain one or any of the substitution strings.

    Args:
        s (str):
            The string into which the replacements are made.
        ip (str):
            The string with which to replace <INSTALL_DIR>.
        dp (str):
            The string with which to replace <DATA_DIR>.
        gp (str):
            The string with which to replace <GLOBAL_DATA>.
        ei (str):
            The string with which to replace <EVENT_ID>.

    Returns:
        str: A new string with the sub-string replacements.
    """

    s = s.replace('<INSTALL_DIR>', ip)
    s = s.replace('<DATA_DIR>', dp)
    s = s.replace('<GLOBAL_DATA>', gp)
    s = s.replace('<EVENT_ID>', ei)
    return s


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
        'nanfloat_type': nanfloat_type,
        'nanfloat_list': nanfloat_list,
        'gmpe_list': gmpe_list,
        'weight_list': weight_list,
        'extent_list': extent_list,
        'status_string': status_string,
    }
    validator = Validator(fdict)
    return validator


def config_error(config, results):
    """
    Parse the results of a ConfigObj validation and log the errors.
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
            logging.error('The "%s" key in the section "%s" failed validation'
                          % (key, ', '.join(section_list)))
            errs += 1
        else:
            logging.error('The following section was missing:%s '
                          % ', '.join(section_list))
            errs += 1
    if errs:
        raise RuntimeError('There %s %d %s in configuration.'
                           % ('was' if errs == 1 else 'were', errs,
                              'error' if errs == 1 else 'errors'))


def check_extra_values(config, logger):
    """
    Checks the config and warns the user if there are any extra entries
    in their config file. This function is based on suggested usage in the
    ConfigObj manual.

    Args:
        config (ConfigObj): A ConfigObj instance.
        logger (logger): The logger to which to write complaints.

    Returns:
        Nothing: Nothing.
    """
    warnings = 0
    for sections, name in get_extra_values(config):

        # this code gets the extra values themselves
        the_section = config
        for section in sections:
            the_section = the_section[section]

        # the_value may be a section or a value
        the_value = the_section[name]

        section_or_value = 'value'
        if isinstance(the_value, dict):
            # Sections are subclasses of dict
            section_or_value = 'section'

        section_string = ', '.join(sections) or "top level"
        logger.warning('Extra entry in section: %s: %s %r is not in spec.' %
                       (section_string, section_or_value, name))
        warnings += 1
    if warnings:
        logger.warning('The extra value(s) may be the result of deprecated '
                       'entries or other changes to the config files; please '
                       'check the conifg files in shakemap/data for the most '
                       'up to date specs.')


def check_config(config, logger):
    """
    Checks that the gmpe, gmice, ipe, and ccf parameters
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


def check_all_configs(configdir):
    data_path = get_data_path()
    specfiles = glob.glob(os.path.join(data_path, '*spec*.conf'))
    missing_files = []
    issues = {}
    exceptions = []
    val = get_custom_validator()
    for tspecfile in specfiles:
        _, specfile = os.path.split(tspecfile)
        configfile = os.path.join(configdir, specfile.replace('spec', ''))
        if not os.path.isfile(configfile):
            missing_files.append(configfile)
            continue
        config = ConfigObj(configfile, configspec=tspecfile,
                           interpolation=False)
        try:
            results = config.validate(val, preserve_errors=True)
        except Exception as e:
            exceptions.append(str(e))
        if not isinstance(results, bool):
            # results = flatten_errors(config, results)
            issues[specfile] = results

    return (missing_files, issues, exceptions)


def nanfloat_type(value):
    """
    Checks to see if value is a float, or NaN, nan, Inf, -Inf, etc.
    Raises a ValidateError exception on failure.

    Args:
        value (str): A string representing a float NaN or Inf.

    Returns:
        float: The input value converted to its float equivalent.

    """
    try:
        out = float(value)
    except ValueError:
        raise ValidateError(value)
    return out


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
        except Exception:
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
            logging.error("list must contain at least %d entr%s" %
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
                logging.error("list must contain at least %d entr%s" %
                              (min, "ies" if int(min) != 1 else "y"))
                raise ValidateError()
        else:
            value = [value]
    if len(value) < int(min):
        logging.error("list must contain at least %d entr%s" %
                      (min, "ies" if int(min) != 1 else "y"))
        raise ValidateError()
    try:
        out = [float(a) for a in value]
    except ValueError:
        logging.error("%s is not a list of floats" % value)
        raise ValidateError()
    np_out = np.array(out)
    if np.any(np_out < 0):
        logging.error("all list values must be >= 0: %s" % value)
        raise ValidateError()
    if len(out) > 0 and np.abs(np.sum(np_out) - 1.0) > 0.01:
        logging.error("weights must sum to 1.0: %s" % value)
        raise ValidateError()

    return out


def nanfloat_list(value, min):
    """
    Checks to see if value is a list of floats, including NaN and Inf.
    Raises a ValidateError exception on failure.

    Args:
        value (str): A string representing a list of floats.

    Returns:
        list: The input string converted to a list of floats.

    """
    min = int(min)
    if isinstance(value, str) and (value == 'None' or value == '[]'):
        value = []
    if isinstance(value, str):
        value = [value]
    if isinstance(value, list) and not value:
        value = []
    if not isinstance(value, list):
        logging.error("'%s' is not a list" % value)
        raise ValidateError()
    if len(value) < min:
        logging.error("extent list must contain %i entries" % min)
        raise ValidateError()
    try:
        out = [float(a) for a in value]
    except ValueError:
        logging.error("%s is not a list of %i floats" % (value, min))
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
        logging.error("'%s' is not a list of at least %s gmpes" % (value, min))
        raise ValidateError()
    for gmpe in value:
        if not isinstance(gmpe, str):
            logging.error("'%s' is not a list of strings" % (value))
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
        logging.error("'%s' is not a list of 4 coordinates" % value)
        raise ValidateError()
    if len(value) != 4:
        logging.error("extent list must contain 4 entries")
        raise ValidateError()
    try:
        out = [float(a) for a in value]
    except ValueError:
        logging.error("%s is not a list of 4 floats" % value)
        raise ValidateError()
    if out[0] < -360.0 or out[0] > 360.0 or \
       out[2] < -360.0 or out[2] > 360.0 or \
       out[1] < -90.0 or out[1] > 90.0 or \
       out[3] < -90.0 or out[3] > 90.0:
        logging.error("Invalid extent: %s "
                      "(-360 <= longitude <= 360, -90 <= latitude <= 90)"
                      % value)
        raise ValidateError()

    return out


def file_type(value):
    """
    Checks to see if value is a valid file or an empty string.
    Raises a ValidateError exception on failure. Does macro substitution
    of <INSTALL_DIR>, <DATA_DIR>, and <GLOBAL_DATA_DIR>.

    Args:
        value (str): A string representing the path to a file.

    Returns:
        str: The input string.

    """
    if not value or value == 'None':
        return ''
    ip, dp = get_config_paths()
    gp = os.path.join(os.path.expanduser('~'), 'shakemap_data')
    value = path_macro_sub(value, ip=ip, dp=dp, gp=gp)
    if not os.path.isfile(value):
        logging.error("file '%s' is not a valid file" % value)
        raise ValidateError(value)
    return value


def directory_type(value):
    """
    Checks to see if value is a valid directory or an empty string.
    Raises a ValidateError exception on failure. Does macro substitution
    of <INSTALL_DIR>, <DATA_DIR>, and <GLOBAL_DATA_DIR>.

    Args:
        value (str): A string representing the path to a directory.

    Returns:
        str: The input string.

    """
    if not value or value == 'None':
        return ''
    ip, dp = get_config_paths()
    gp = os.path.join(os.path.expanduser('~'), 'shakemap_data')
    value = path_macro_sub(value, ip=ip, dp=dp, gp=gp)
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
        logging.error("'%s' is not a list of at least 1 float" % (value))
        raise ValidateError()
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, list) or len(value) < 1:
        logging.error("'%s' is not a list of at least 1 float" % (value))
        raise ValidateError()
    fvalue = []
    for val in value:
        try:
            fval = float(val)
        except ValueError:
            logging.error("'%s' is not a list of floats" % (value))
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
        logging.error("'%s' is not a float" % (value))
        raise ValidateError()
    try:
        fval = float(value)
    except ValueError:
        logging.error("'%s' is not a float" % (value))
        raise ValidateError()
    return fval


def cfg_bool(value):
    """
    Converts (if possible) the input string to a bool. Raises
    ValidateError if the input can't be converted to a bool.

    Args:
        value (str): A string to be converted to a bool.

    Returns:
        bool: The input converted to a bool.

    Raises:
        ValidateError
    """
    if not isinstance(value, (str, bool)) or not value or value == 'None':
        logging.error("'%s' is not a bool" % (value))
        raise ValidateError()
    try:
        if value.lower() in ['true', 't', 'yes', 'y', '1']:
            bval = True
        else:
            bval = False
    except ValueError:
        logging.error("'%s' is not a bool" % (value))
        raise ValidateError()
    return bval


def check_profile_config(config):
    """
    Validation checks on the profile config. At least one profile must exist
    (otherwise exit) and the paths for each profile should exist, otherwise the
    profile entry is removed.

    Args:
        config (ConfigObj): The ConfigObj instance.
    """
    # Check that at least one profile exists
    if 'profiles' not in config:
        logging.error('There are currently no profiles. Use "sm_profile '
                      '-c <profile>" to create one.')
        sys.exit(1)
    # Check that the paths for each profile exist
    for profile in config['profiles'].keys():
        data_exists = os.path.isdir(config['profiles'][profile]['data_path'])
        delete_profile = False
        if not data_exists:
            logging.warn('Data path for profile %s does not exist.' % profile)
            delete_profile = True
        install_exists = os.path.isdir(
            config['profiles'][profile]['install_path'])
        if not install_exists:
            logging.warn(
                'Install path for profile %s does not exist.' % profile)
            delete_profile = True
        if delete_profile:
            logging.warn('    Deleting profile %s.' % profile)
            del config['profiles'][profile]
            config.write()
    return config
