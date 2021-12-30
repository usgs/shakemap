import sys
import os
from importlib import import_module
from urllib import request
from urllib.error import HTTPError
import json
import logging
from collections import OrderedDict
from shakemap.utils.config import get_config_paths
import time

# third party imports
from configobj import ConfigObj

# url template for json file describing everything we know about a network
NETWORK_TEMPLATE = (
    "https://earthquake.usgs.gov/data/comcat/" "contributor/[NETID]/index.json"
)


def migrate_gmpe(old_gmpe, config=None):
    """Return the GMPE that should be used to replace SM3.5 GMPE

    By default, this uses the migrate.conf file found in the ShakeMap
    repository. Users can optionally pass in their own config.

    Args:
        old_gmpe (str):
            ShakeMap 3.5 GMPE string
        config (dict):
            Input configobj dict or None.

    Returns:
        tuple: New GMPE string, and GMPE reference string.
    """
    if config is None:
        install_path, _ = get_config_paths()
        if not os.path.isdir(install_path):
            raise OSError(f"{install_path} is not a valid directory.")
        config_file = os.path.join(install_path, "config", "migrate.conf")
        if os.path.isfile(config_file):
            config = ConfigObj(config_file)
        else:
            raise OSError(f"{config_file} not found.")
    if old_gmpe not in config["modules"]:
        raise KeyError(f"ShakeMap 3.5 GMPE {old_gmpe} not found in migrate.conf.")
    new_gmpe = config["modules"][old_gmpe]["openquake"]
    reference = config["modules"][old_gmpe]["reference"]
    return (new_gmpe, reference)


def set_gmpe(gmpe, config, eventid):
    gmpe_list = [gmpe]
    weight_list = [1.0]
    gmpe_set = "gmpe_" + eventid + "_custom"
    config["modeling"]["gmpe"] = gmpe_set
    config["gmpe_sets"] = OrderedDict(
        [
            (
                gmpe_set,
                OrderedDict(
                    [
                        ("gmpes", gmpe_list),
                        ("weights", weight_list),
                        ("weights_large_dist", "None"),
                        ("dist_cutoff", "nan"),
                        ("site_gmpes", "None"),
                        ("weights_site_gmpes", "None"),
                    ]
                ),
            )
        ]
    )
    return config


def get_network_name(netid):
    """Return a string representing a name of a network given its ID.

    Note: Uses an internet connection to ComCat.

    Args:
        netid (str): Usually two-character network ID (us,ci, etc.)
    Returns:
        str: Network name, or "unknown" if input netid is invalid.
    """
    url = NETWORK_TEMPLATE.replace("[NETID]", netid)
    network = "unknown"
    fails = 0
    while fails < 3:
        try:
            fh = request.urlopen(url)
            data = fh.read().decode("utf-8")
            jdict = json.loads(data)
            fh.close()
            network = jdict["title"]
        except HTTPError:
            error_str = (
                "No network description found for %s. You may "
                "want to make sure that the pages of the form %s "
                "still exist and contact the ShakeMap developers "
                "if they do not." % (netid, NETWORK_TEMPLATE)
            )
            logging.warning(error_str)
            fails += 1
        except Exception as e:
            logging.warning(f"Error in get_network_name: {str(e)}")
            logging.warning("Will try %d more times" % (3 - fails))
            fails += 1
            time.sleep(20)
        else:
            break

    return network


def get_object_from_config(obj, section, cfg, *args):
    """
    Helper function for things (ipe, gmice, ccf) that don't have a
    fromConfig() constructor yet. Instantiates an instance of a
    class from config entry, 'name', that has a corresponding
    'name'_module dictionary of class name, module path.

    Args:
        obj (str):
            Name of the parameter in the config file that specifies
            the object to be instantiated.
        section (str):
            The section of the config in which 'obj' resides.
        cfg (dict):
            The configuration file in which 'obj' and the
            module definitions reside.
        args:
            Additional arguments that will be passed to the constructor
            of the thing being instantiated.

    Returns:
        object: An instance of the object specified by the 'obj' parameter
            in the config file..
    """
    cls_abbr = cfg[section][obj]
    mods = obj + "_modules"
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
        bool: The "answer" return value is True for "yes" or False for "no".

    """
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError(f"invalid default answer: '{default}'")

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")
