#!/usr/bin/env python

import subprocess


def get_command_output(cmd):
    """
    Method for calling external system command.

    :param cmd: 
        String command (e.g., 'ls -l', etc.). 
    :returns: 
        Three-element tuple containing a boolean indicating success or failure, 
        the stdout from running the command, and stderr.
    """
    proc = subprocess.Popen(cmd,
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE
                            )
    stdout, stderr = proc.communicate()
    retcode = proc.returncode
    if retcode == 0:
        retcode = True
    else:
        retcode = False
    return (retcode, stdout, stderr)
