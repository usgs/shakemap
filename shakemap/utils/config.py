import os.path
import pkg_resources

# third party libraries
import numpy as np
from configobj import ConfigObj, flatten_errors
from validate import Validator, ValidateError

def get_data_path():
    """
    Return the path to the shakemap package data directory holding
    configspec.conf, the template config files, and the example config
    files.

    Returns:
        (str): The full path to the data directory.

    """
    return pkg_resources.resource_filename('shakemap', 'data')

def get_configspec():
    """
    Returns the full path to the ShakeMap configspec.conf file.

    Returns:
        (str): The path to the configspec.

    """
    return os.path.join(get_data_path(), 'configspec.conf')


def get_config_paths():
    """
    Returns two paths based on the currently selected profile in the 
    user's ~/.shakemap/profile.conf: 1) the path to the ShakeMap
    installation directory; 2) the path to the data directory.

    Returns:
        (str, str): The paths to the ShakeMap install directory
        and the data directory.
    """
    config_file = os.path.join(os.path.expanduser('~'), '.shakemap', 
                               'profiles.conf')
    config = ConfigObj(config_file)
    profile_name = config['profile']
    profile = config['profiles'][profile_name]
    install = profile['install_path']
    data = profile['data_path']
    return (install,data)

def get_custom_validator():
    """
    Returns a validator suitable for validating the ShakeMap config
    files.

    Returns:
        (:class:`Validator`): A Validator object.

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
        (Nothing): Nothing

    """
    errs = 0
    for (section_list, key, _) in flatten_errors(config, results):
        if key is not None:
            print('The "%s" key in the section "%s" failed validation' % 
                    (key, ', '.join(section_list)))
            errs += 1
        else:
            print('The following section was missing:%s ' % 
                    ', '.join(section_list))
            errs += 1
    if errs:
        raise RuntimeError('There %s %d %s in configuration.' %
                    ('was' if errs == 1 else 'were', errs, 
                     'error' if errs == 1 else 'errors'))

def check_config(config):
    """
    Checks that the gmpe, gmice, ipe, ccf, and component parameters
    in config are defined in their respective sections. Raises a 
    ValidateError exception if an error is encountered.

    Args:
        config (ConfigObj): A ConfigObj instance.

    Returns:
        (Nothing): Nothing.

    """
    if config['modeling']['gmpe'] not in config['gmpe_sets']:
        print('Configuration error: gmpe %s not in gmpe_sets' %
              (config['modeling']['gmpe']))
        raise ValidateError()
    if config['modeling']['gmice'] not in config['gmice_modules']:
        print('Configuration error: gmice %s not in gmice_modules' %
              (config['modeling']['gmice']))
        raise ValidateError()
    if config['modeling']['ipe'] not in config['ipe_modules']:
        print('Configuration error: ipe %s not in ipe_modules' %
              (config['modeling']['ipe']))
        raise ValidateError()
    if config['modeling']['ccf'] not in config['ccf_modules']:
        print('Configuration error: ccf %s not in ccf_modules' %
              (config['modeling']['ccf']))
        raise ValidateError()
    if config['interp']['component'] not in config['component_modules']:
        print('Configuration error: component %s not in component_modules' %
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
        (float): The input value converted to decimal degrees.

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
        (list): The input string converted to a list of floats.

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
    if len(out) > 0 and np.sum(np_out) != 1:
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
        (list): The input string converted to a list of GMPEs.

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
        (list): The input string converted to a list of floats.

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
        (string): The input string.

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
        (string): The input string.

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
        (string): The input string. 'automatic' is returned if
        value is empty.

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
        (list): The input converted to a list of floats.

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
        (float): The input converted to a float.

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
