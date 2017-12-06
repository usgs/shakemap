import sys
from importlib import import_module


def get_object_from_config(obj, section, cfg, *args):
    """
    Helper function for things (ipe, gmice, ccf) that don't have a
    fromConfig() constructor yet. Instantiates an instance of a
    class from config entry, 'name', that has a corresponding
    'name'_module dictionary of class name, module path.

    Args:
        obj (str): Name of the parameter in the config file that specifies
            the object to be instantiated.
        section (str): The section of the config in which 'obj' resides.
        cfg (dictionary): The configuration file in which 'obj' and the
            module definitions reside.
        args: Additional arguments that will be passed to the constructor
            of the thing being instantiated.

    Returns:
        object: An instance of the object specified by the 'obj' parameter
            in the config file..
    """
    cls_abbr = cfg[section][obj]
    mods = obj + '_modules'
    (cname, mpath) = cfg[mods][cls_abbr]
    return getattr(import_module(mpath), cname)(*args)


#
# This is from stackoverflow where it was reproduced from Recipe 577058
#
def query_yes_no(question, default="yes"):
    """
    Ask a yes/no question via raw_input() and return their answer.

    Args:
        question (str): a string that is presented to the user.
        default (str): the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).

    Returns:
        (bool): The "answer" return value is True for "yes" or False for "no".

    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def path_macro_sub(s, ip, dp):
    """
    Replace macros with current paths:

    - <INSTALL_DIR> is replaced with the contents of ip
    - <DATA_DIR> is replaced with the contents of dp

    E.g., path_macro_sub("<INSTALL_DIR>/<DATA_DIR>", "hello", "world") 
    would return "hello/world". It is not an error if the original string 
    does not contain one or both of <INSTALL_DIR> or <DATA_DIR>.

    Args:
        s (str): The string into which the replacements are made.
        ip (str): The string with which to replace <INSTALL_DIR>.
        dp (str): The string with which to replace <DATA_DIR>.

    Returns:
        (str): A new string with the sub-string replacements.
    """

    return s.replace('<INSTALL_DIR>', ip).replace('<DATA_DIR>', dp)
