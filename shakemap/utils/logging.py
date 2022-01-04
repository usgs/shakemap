# stdlib imports
import os
import os.path
import logging
import logging.config
import logging.handlers

# third party libraries
from configobj import ConfigObj
from validate import Validator

# local libraries
from shakemap.utils.config import get_config_paths, get_configspec, config_error

REQ_FIELDS = {
    "logging.handlers.TimedRotatingFileHandler": [
        "level",
        "formatter",
        "class",
        "when",
        "filename",
    ],
    "logging.FileHandler": ["level", "formatter", "class", "filename"],
    "logging.handlers.SMTPHandler": [
        "level",
        "formatter",
        "mailhost",
        "fromaddr",
        "toaddrs",
        "subject",
        "class",
    ],
}

# mapping of string versions of logging module logging levels
# to the corresponding constants.
LOG_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
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
        fmt = config["formatters"]["standard"]["format"]
        datefmt = config["formatters"]["standard"]["datefmt"]
        # create a console handler, with verbosity setting chosen by user
        if log_option == "debug":
            level = logging.DEBUG
        elif log_option == "quiet":
            level = logging.ERROR
        else:  # default interactive
            level = logging.INFO

        logdict = {
            "version": 1,
            "formatters": {"standard": {"format": fmt, "datefmt": datefmt}},
            "handlers": {
                "stream": {
                    "level": level,
                    "formatter": "standard",
                    "class": "logging.StreamHandler",
                }
            },
            "loggers": {
                "": {"handlers": ["stream"], "level": level, "propagate": True}
            },
        }

        logging.config.dictConfig(logdict)
    else:
        event_log_dir = os.path.join(data_path, eventid)
        if not os.path.isdir(event_log_dir):
            raise NotADirectoryError(f"Can't open log file: event {eventid} not found")
        event_log_file = os.path.join(event_log_dir, "shake.log")
        config["handlers"]["event_file"]["filename"] = event_log_file

        if log_option == "debug":
            config["handlers"]["event_file"]["level"] = logging.DEBUG
        elif log_option == "quiet":
            config["handlers"]["event_file"]["level"] = logging.ERROR
        elif (
            "level" in config["handlers"]["event_file"]
            and config["handlers"]["event_file"]["level"] in LOG_LEVELS
        ):
            config["handlers"]["event_file"]["level"] = LOG_LEVELS[
                config["handlers"]["event_file"]["level"]
            ]
        else:
            config["handlers"]["event_file"]["level"] = logging.INFO
        global_log_dir = os.path.join(install_path, "logs")

        if not os.path.isdir(global_log_dir):
            os.makedirs(global_log_dir, exist_ok=True)
        global_log_file = os.path.join(global_log_dir, "shake.log")
        config["handlers"]["global_file"]["filename"] = global_log_file

        if log_option == "debug":
            config["handlers"]["global_file"]["level"] = logging.DEBUG
        elif log_option == "quiet":
            config["handlers"]["global_file"]["level"] = logging.ERROR
        elif (
            "level" in config["handlers"]["global_file"]
            and config["handlers"]["global_file"]["level"] in LOG_LEVELS
        ):
            config["handlers"]["global_file"]["level"] = LOG_LEVELS[
                config["handlers"]["global_file"]["level"]
            ]
        else:
            config["handlers"]["global_file"]["level"] = logging.INFO

        if log_option == "debug":
            config["loggers"][""]["level"] = logging.DEBUG
        elif log_option == "quiet":
            config["loggers"][""]["level"] = logging.ERROR
        elif (
            "level" in config["loggers"][""]
            and config["loggers"][""]["level"] in LOG_LEVELS
        ):
            config["loggers"][""]["level"] = LOG_LEVELS[config["loggers"][""]["level"]]
        else:
            config["loggers"][""]["level"] = logging.INFO

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

    install_path, _ = get_config_paths()
    conf_file = os.path.join(install_path, "config", "logging.conf")
    spec_file = get_configspec(config="logging")
    log_config = ConfigObj(conf_file, configspec=spec_file, interpolation="template")

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
    log_name = log_config["loggers"].keys()[0]
    log_config["loggers"][""] = log_config["loggers"][log_name]
    del log_config["loggers"][log_name]
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
    for handler in config["handlers"].values():
        myclass = handler["class"]
        req_fields = REQ_FIELDS[myclass]
        for key in handler.keys():
            if key not in req_fields:
                del handler[key]


def get_generic_logger(logfile=None, fmt=None, datefmt=None, level=None):
    """Returns a generic logging configuration dictionary that may be used
    by programs like receive_amps and receive_origins to create a logger.
    If the fmt and/or the  datefmt are not specified, they will be taken
    from the config in the operator's logging.conf. If logfile is not
    specified, the configuration will be for a stream handler. If level
    is not specified, it will be logging.INFO.

    Args:
        logfile:
            The path to the file to receive the logging information. This
            will be a TimedRotatingFilehandler that will reset at midnight
            every night. If logfile is not specified, a StreamHandler will
            be configured.
        fmt:
            The format string for the logging messages. If it is not
            supplied, the string from the current profile's logging.conf
            will be used.
        datefmt:
            The format string for the logging date/time. If it is not
            supplied, the string from the current profile's logging.conf
            will be used.
        level:
            The logging level. If not specified, it will be 'INFO'. Other
            valid choices are 'DEBUG', 'WARNING', or 'ERROR'.

    Returns:
        logger: A logger suitable for logging messages.
    """
    if fmt is None or datefmt is None:
        config = get_logging_config()
    if fmt is None:
        fmt = config["formatters"]["standard"]["format"]
    if datefmt is None:
        datefmt = config["formatters"]["standard"]["datefmt"]

    if logfile is None:
        handler = "logging.StreamHandler"
    else:
        handler = "logging.handlers.TimedRotatingFileHandler"

    if level is None:
        level = "INFO"

    log_cfg = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "standard": {
                "format": fmt,
                "datefmt": datefmt,
            },
        },
        "handlers": {
            "default": {
                "level": "INFO",
                "formatter": "standard",
                "class": handler,
                "when": "midnight",
            },
        },
        "loggers": {
            "": {"handlers": ["default"], "level": level, "propagate": True},
        },
    }

    if logfile:
        log_cfg["handlers"]["default"]["filename"] = logfile
    logging.config.dictConfig(log_cfg)
    logger = logging.getLogger()
    return logger
