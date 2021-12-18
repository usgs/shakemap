#!/usr/bin/env python

import os
import os.path
import glob

# import pytest
import importlib.util
from tempfile import mkstemp
import shutil
from pathlib import Path

from shakemap.utils.dependencies import sha_sum, DepNode, CommandDatabase

homedir = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.abspath(os.path.join(homedir, "data", "dependencies"))


def replace(file_path, pattern, subst):
    """Function I found on stackoverflow to replace a pattern in a file"""
    # Create temp file
    fh, abs_path = mkstemp()
    with os.fdopen(fh, "w") as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    # Remove original file
    os.remove(file_path)
    # Move new file
    shutil.move(abs_path, file_path)


def test_sha256():
    filename = os.path.join(datadir, "data", "sha256_test.txt")
    assert (
        sha_sum(filename)
        == "7309dd56bf7e1c980482350aeed5c04040315fb1ece0dba25673042edcb1eb25"
    )


def test_depnode():

    # Build a dependency tree
    cmd0 = "Command0"
    cmd1 = "Command1"
    cmd2 = "Command2"
    cmd3 = "Command3"
    cmd4 = "Command4"
    cmd5 = "Command5"
    cmd6 = "Command6"
    cmd7 = "Command7"
    deps0 = [("file1", 0), ("file2", 1)]

    root = DepNode(deps0, cmd0)  # level 0

    leaves = root.getLeaves()
    assert len(leaves) == 1
    assert root in leaves

    parent1 = DepNode(deps0, cmd1, root)  # level 1
    parent2 = DepNode(deps0, cmd2, root)  # level 1
    parent11 = DepNode(deps0, cmd3, parent1)  # level 2
    parent12 = DepNode(deps0, cmd4, parent1)  # level 2
    parent21 = DepNode(deps0, cmd5, parent2)  # level 2
    parent22 = DepNode(deps0, cmd6, parent2)  # level 2
    parent211 = DepNode(deps0, cmd7, parent21)  # level 3

    # Check that the final leaf has reasonable values
    assert parent211.parents() == []
    assert parent211.child() == parent21
    assert parent211.deps() == deps0
    assert parent211.cmd() == cmd7
    assert parent211.level() == 3
    assert parent211.child().child().cmd() == cmd2

    # Check that the leaves (nodes with no parents) are as expected
    leaves = root.getLeaves()
    assert len(leaves) == 4
    assert parent11 in leaves
    assert parent12 in leaves
    assert parent22 in leaves
    assert parent211 in leaves

    # Check that some of the parents are properly organized
    parents = parent1.parents()
    assert len(parents) == 2
    assert parent11 in parents
    assert parent12 in parents

    parents = root.parents()
    assert len(parents) == 2
    assert parent1 in parents
    assert parent2 in parents

    # Check that the levels of all the nodes are what they should be
    assert parent1.level() == 1
    assert parent2.level() == 1
    assert parent11.level() == 2
    assert parent12.level() == 2
    assert parent21.level() == 2
    assert parent22.level() == 2


file_list = {
    "config": ["module1a.conf", "module2a.conf", "module1b.conf", "module2b.conf"],
    os.path.join("dummy_event", "current"): ["module1_dep1.txt", "module1_dep2.txt"],
    "modules": ["module1.py", "module2.py"],
}


def make_files(datadir):
    for subdir, filelist in file_list.items():
        dest = os.path.join(datadir, subdir)
        src = os.path.join(dest, "src")
        for fp in filelist:
            shutil.copy(os.path.join(src, fp), os.path.join(dest, fp))


def clear_files(datadir):
    for subdir, filelist in file_list.items():
        dest = os.path.join(datadir, subdir)
        for fp in filelist:
            if os.path.isfile(os.path.join(dest, fp)):
                os.remove(os.path.join(dest, fp))
    # Also get rid of the targets:
    t1 = os.path.join(datadir, "dummy_event", "current", "module1_target.txt")
    if os.path.isfile(t1):
        os.remove(t1)
    t2 = os.path.join(
        datadir, "dummy_event", "current", "products", "module2_target.txt"
    )
    if os.path.isfile(t2):
        os.remove(t2)


def test_command_database():

    clear_files(datadir)
    make_files(datadir)

    dbfiles = glob.glob(os.path.join(datadir, "data", "dependencies.*"))
    for dbfile in dbfiles:
        os.remove(dbfile)

    eventdir = os.path.join(datadir, "dummy_event", "current")
    mod1_target = os.path.join(eventdir, "module1_target.txt")
    mod2_target = os.path.join(eventdir, "products", "module2_target.txt")

    # Import the modules
    moddir = os.path.join(datadir, "modules")
    # Module 1
    mod1_path = os.path.join(moddir, "module1.py")
    spec = importlib.util.spec_from_file_location("module1", mod1_path)
    mod1 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod1)
    module1 = mod1.Module1()

    # Module 2
    mod2_path = os.path.join(moddir, "module2.py")
    spec = importlib.util.spec_from_file_location("module2", mod2_path)
    mod2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod2)
    module2 = mod2.Module2()

    classes = {
        module1.command_name: {"class": module1, "mfile": mod1_path},
        module2.command_name: {"class": module2, "mfile": mod2_path},
    }

    #
    # This is the test setup:
    # command: module1
    #   dependencies: module1_dep1.txt, module1_dep2.txt
    #   configs: module1a.conf, module1b.conf
    #   target: module1_target.txt
    # command: module2
    #   dependencies: module1_target.txt
    #   configs: module2a.conf, module2b.conf
    #   target: module2_target.txt
    #
    cdb = CommandDatabase(classes, "dummy_event", datadir, datadir)
    #
    # We start with everything out of date, and then bring everything
    # up to date and step through the ways things can change
    #
    status = cdb.getDependencyStatus("module2")
    assert status == ["module1", "module2"]
    cdb.updateCommand("module1")

    # Ran module1, but it's target isn't there yet
    status = cdb.getDependencyStatus("module2")
    assert status == ["module1", "module2"]

    # Now make the target, but it's out of date vis a vis module2
    Path(mod1_target).touch()
    status = cdb.getDependencyStatus("module2")
    assert status == ["module2"]

    # Ran module2, but it's target isn't there yet
    cdb.updateCommand("module2")
    status = cdb.getDependencyStatus("module2")
    assert status == ["module2"]

    # Now make the target for module2
    Path(mod2_target).touch()
    status = cdb.getDependencyStatus("module2")
    assert status == []

    # Change a moodule1 config file:
    configdir = os.path.join(datadir, "config")
    mod1_config = os.path.join(configdir, "module1a.conf")
    replace(mod1_config, "# Comment1", "# Comment2")
    status = cdb.getDependencyStatus("module2")
    assert status == ["module1", "module2"]

    # Update the command, but since the target didn't change, module2 is
    # also up to date
    cdb.updateCommand("module1")
    status = cdb.getDependencyStatus("module2")
    assert status == []

    # Change a moodule2 config file:
    mod2_config = os.path.join(configdir, "module2a.conf")
    replace(mod2_config, "# Comment1", "# Comment2")
    status = cdb.getDependencyStatus("module2")
    assert status == ["module2"]

    # Change it back
    replace(mod2_config, "# Comment2", "# Comment1")
    status = cdb.getDependencyStatus("module2")
    assert status == []

    # Close and reopen the database: everything should be okay
    cdb.close()
    cdb = CommandDatabase(classes, "dummy_event", datadir, datadir)
    status = cdb.getDependencyStatus("module2")
    assert status == []

    # Change module1 source; this requires a close and reopen of
    # the database since the command checksums are done at startup
    moddir = os.path.join(datadir, "modules")
    mod1_src = os.path.join(moddir, "module1.py")
    replace(mod1_src, "# Comment1", "# Comment2")
    cdb.close()
    cdb = CommandDatabase(classes, "dummy_event", datadir, datadir)
    status = cdb.getDependencyStatus("module2")
    assert status == ["module1", "module2"]

    # Update the command. Again, since the target didn't change,
    # module2 is up to date.
    cdb.updateCommand("module1")
    status = cdb.getDependencyStatus("module2")
    assert status == []

    # Remove a required dependency; maybe this should throw an exception?
    os.remove(os.path.join(eventdir, "module1_dep1.txt"))
    os.remove(os.path.join(eventdir, "module1_dep2.txt"))
    status = cdb.getDependencyStatus("module2")
    assert status == ["module1", "module2"]

    dbfiles = glob.glob(os.path.join(datadir, "data", "dependencies.*"))
    for dbfile in dbfiles:
        os.remove(dbfile)

    clear_files(datadir)


if __name__ == "__main__":
    test_sha256()
    test_depnode()
    test_command_database()
