import os.path
import hashlib
import sqlite3
import glob
import re
import logging


from shakemap.utils.config import get_config_paths


def sha_sum(fname):
    """Returns the sha256 checksum of a file.

    Args:
        fname (str): Path to a file.

    Returns:
        str: The sha256 checksum as a hex string.
    """
    hasher = hashlib.sha256()
    with open(fname, "rb") as fd:
        for block in iter(lambda: fd.read(65536), b""):
            hasher.update(block)
    return hasher.hexdigest()


class DepNode(object):
    """A helper class to create a dependency tree made up of "nodes" with
    parents, children, and data.
    """

    def __init__(self, deps, cmd, child=None):
        """Create a new node. When called without a child node, this will be
        the root node of a tree. When called with a child, this node will
        be a parent of the child node.

        Args:
            deps (list): A list of (file, required (bool)) that are the
                dependencies of cmd.
            cmd (str): The name of the ShakeMap coremod command.
            child (class DepNode): The child of the created node; if None,
                the node will be the root node of a tree.
        """
        self._deps = deps
        self._cmd = cmd
        self._child = child
        self._parents = []
        self._level = 0
        if child is not None:
            self._level = child._level + 1
            child._parents.append(self)

    def parents(self):
        """Returns the parents of the calling node.

        Args:
            none

        Returns:
            list: Returns a list of parent nodes of the calling node.
        """
        return self._parents

    def child(self):
        """Returns the child of the calling node.

        Args:
            none

        Returns:
            class DepNode or None: Returns the child node of the calling node
            or None if the calling node is the root.
        """
        return self._child

    def deps(self):
        """Returns the dependency list of the calling node.

        Args:
            none

        Returns:
            list: A list of (file (str), required (bool)) tuples that are the
                dependencies of the node's cmd.
        """
        return self._deps

    def cmd(self):
        """Returns the command stored in the calling node.

        Args:
            none

        Returns:
            str: The command stored in the node.
        """
        return self._cmd

    def level(self):
        """Returns the tree level of the calling node.

        Args:
            none

        Returns:
            int: The tree level of the node.
        """
        return self._level

    def getLeaves(self):
        """Returns all the leaf nodes (i.e., nodes with no parents) in the
        tree that starts at the calling node.

        Args:
            none

        Returns:
            list: A list of leaf nodes (i.e., nodes with no parents) in the
            tree that starts at the calling node.
        """
        if len(self._parents) == 0:
            return [self]
        leaves = []
        for parent in self._parents:
            if len(parent._parents) == 0:
                leaves += [parent]
                continue
            leaves += parent.getLeaves()
        return leaves


class CommandDatabase(object):
    """A class to keep track of dependencies."""

    def __init__(self, classes, eventid, install_path=None, data_path=None):
        """
        Create a new (or open an exsiting) dependency database for checking
        the status of command dependencies.

        Args:
            classes (dict): A dictionary commands and their associated data
                as produced by the 'shake' program's _get_command_classes()
                function.
            eventid (str): The event ID of the event being processed.

        Returns:
            class CommandDatabase: An instance of the CommandDatabase class.
        """
        if install_path is None or data_path is None:
            install_path, data_path = get_config_paths()
        event_path = os.path.join(data_path, eventid, "current")

        self.eventid = eventid
        self.event_path = event_path
        #
        # The "c" database is a temporary in-memory database just for
        # this run of a particular event
        #
        self.cconnect = con = sqlite3.connect(":memory:")
        self.ccursor = cur = con.cursor()

        cur.execute(
            """CREATE TABLE dependencies
                       (command text, dependency text, required integer)"""
        )
        con.commit()

        config_path = os.path.join(install_path, "config")
        self.targets = {}
        self.cmdsums = {}
        dependencies = []
        for cmd, cd in classes.items():
            if cd["class"].targets is not None:
                for target in cd["class"].targets:
                    self.targets[cmd] = self.targets.get(cmd, []) + [
                        re.compile(os.path.join(event_path, target))
                    ]
            else:
                self.targets[cmd] = None
            if cd["class"].dependencies is not None:
                for depend in cd["class"].dependencies:
                    dependencies.append(
                        (cmd, os.path.join(event_path, depend[0]), int(depend[1]))
                    )
            # Configs are dependencies, too, they just live in a different
            # directory
            if cd["class"].configs is not None:
                for config in cd["class"].configs:
                    dependencies.append((cmd, os.path.join(config_path, config), 1))

            self.cmdsums[cmd] = sha_sum(cd["mfile"])

        if len(dependencies) > 0:
            cur.executemany("INSERT into dependencies values (?, ?, ?)", dependencies)
            con.commit()

        #
        # The dependencies.db database keeps track of the state of the
        # dependencies from previous runs
        #
        dbfile = os.path.join(install_path, "data", "dependencies.db")
        if not os.path.isfile(dbfile):
            # It doesn't exist, so create it
            # TODO: Better checking would probably be good here
            self.fconnect = sqlite3.connect(dbfile, timeout=15)
            self.fcursor = self.fconnect.cursor()
            self.fcursor.execute(
                """CREATE TABLE file_checksums (command text NOT NULL,
                   file text NOT NULL,
                   checksum text, PRIMARY KEY (command, file))"""
            )
            self.fcursor.execute(
                """CREATE TABLE cmd_checksums (command text NOT NULL,
                   eventid text NOT NULL,
                   checksum text, PRIMARY KEY (command, eventid))"""
            )
            self.fconnect.commit()
        else:
            self.fconnect = sqlite3.connect(dbfile, timeout=15)
            self.fcursor = self.fconnect.cursor()
        self.fcursor.execute("PRAGMA journal_mode = WAL")

    def close(self):
        """Closes the object and associated databases"""
        self.fcursor.close()
        self.fconnect.close()
        self.ccursor.close()
        self.cconnect.close()

    def getDependencyStatus(self, cmd):
        """Return the status of the specified command and its ancestors.

        Args:
            cmd (str): The command name of one of the configured modules
                for shake (i.e., the "coremods")

        Returns:
            list: A list of tuples of type (command, status) where the
            status is one of:

                missing-depencencies: the command cannot run because of
                    missing dependencies
                out-of-date: the dependencies are newer than the output, so
                    the command should be run
                okay: The output is up to date relative to the dependencies

        """
        tree = self.__buildDependencyTree(cmd)

        status = self.__getDependencyStatus(tree)

        return status

    def __buildDependencyTree(self, cmd, root=None):
        """Internal function to build the dependency tree for a command."""
        deps = self.__getDependencies(cmd)
        child = root = DepNode(deps, cmd, root)
        for fp in deps:
            fp_exists = len(glob.glob(fp[0]))
            # If a file is optional and does not exist, we don't want to
            # force the system to make it. (E.g., we don't want to force
            # 'select' to run if model_select.conf isn't already there --
            # that's the user's choice.
            if fp[1] == 0 and not fp_exists:
                continue
            parent_cmd = self.__findCmd(fp[0])
            if parent_cmd:
                self.__buildDependencyTree(parent_cmd, child)
        return root

    def __getDependencyStatus(self, tree):
        """Internal function to create a list of out-of-date commands."""
        status = []
        for leaf in tree.getLeaves():
            status += self.__traceLeaf(leaf)

        # Make a list of unique commands, in order from ancestors to
        # descendents: set() gets the unique elements, sorted() puts
        # them in proper order based on tree "level", and the list
        # comprehension extracts the command name into a list.
        status = [x[0] for x in sorted(set(status), key=lambda x: x[1], reverse=True)]
        return status

    def __traceLeaf(self, leaf):
        """Internal function to produce a list of all commands that are
        descendents of the leaf argument's command.
        """
        clean = True
        status = []
        while leaf:
            #
            # If any command is out of date, all of its descendents are also
            # out of date -- no need to check their status
            #
            if len(status) > 0:
                status.append((leaf.cmd(), leaf.level()))
                leaf = leaf.child()
                continue
            #
            # If command is out of date, all descendent commands need to be
            # rerun; otherwise ok
            #
            sql = "SELECT checksum FROM cmd_checksums WHERE command=? AND " "eventid=?"
            self.fcursor.execute(sql, (leaf.cmd(), self.eventid))
            result = self.fcursor.fetchall()

            if len(result) == 0 or result[0][0] != self.cmdsums[leaf.cmd()]:
                # command isn't in DB, or file checksum has
                # changed: need to rerun
                status.append((leaf.cmd(), leaf.level()))
                leaf = leaf.child()
                continue
            #
            # If target is missing, the command needs to be run
            # This is only semi-reliable because it only checks that
            # at least one of each target pattern is there.
            #
            explist = self.targets[leaf.cmd()]
            if explist is not None:
                for exp in explist:
                    found_target = False
                    files = glob.glob(os.path.join(self.event_path, "*"))
                    files += glob.glob(os.path.join(self.event_path, "products", "*"))
                    for filename in files:
                        if exp.fullmatch(filename):
                            found_target = True
                            break
                    if found_target is False:
                        status.append((leaf.cmd(), leaf.level()))
                        leaf = leaf.child()
                        clean = False
                        break
            if clean is False:
                continue

            #
            # If dep file is out of date, child command and all descendent
            # commands need to be rerun
            #
            sql = "SELECT checksum FROM file_checksums WHERE command=? " "AND file=?"
            for pattern, required in leaf.deps():
                found = glob.glob(pattern)
                if len(found) == 0 and required:
                    # Didn't find a required file -- this prevents a run
                    logging.warning(f"Did not find required dependency: {pattern}")
                    clean = False
                    break
                if len(found) > 0:
                    for fp in found:
                        self.fcursor.execute(sql, (leaf.cmd(), fp))
                        result = self.fcursor.fetchall()
                        if len(result) == 0 or result[0][0] != sha_sum(fp):
                            # File is a dependency, but isn't in DB, or
                            # file checksum has changed: need to rerun
                            clean = False
                            break
                if clean is False:
                    break

            if clean is False:
                status.append((leaf.cmd(), leaf.level()))

            leaf = leaf.child()

        return status

    def updateCommand(self, cmd):
        """Update the command and dependency checksums of a command that
        has just been run.

        Args:
            cmd (str): The name of the command that has just been run.

        Returns:
            nothing: Nothing.
        """
        sql = (
            "INSERT OR REPLACE INTO file_checksums (command, file, "
            "checksum) values (?, ?, ?)"
        )
        args = []
        files = self.__getDependencies(cmd)
        for fp in files:
            # Only do dependencies that actually exist
            efiles = glob.glob(fp[0])
            for ef in efiles:
                csum = sha_sum(ef)
                args.append((cmd, ef, csum))
        if len(args) > 0:
            self.fcursor.executemany(sql, args)

        sql = (
            "INSERT OR REPLACE INTO cmd_checksums (command, eventid, "
            "checksum) values (?, ?, ?)"
        )
        self.fcursor.execute(sql, (cmd, self.eventid, self.cmdsums[cmd]))
        self.fconnect.commit()
        return

    def __getDependencies(self, cmd):
        """Internal function to get the file dependencies of a command."""
        self.ccursor.execute(
            "SELECT dependency, required FROM dependencies WHERE command=?", (cmd,)
        )
        return self.ccursor.fetchall()

    def __findCmd(self, filepath):
        """Find the command that produces the target file named by
        filepath.
        """
        for cmd, explist in self.targets.items():
            if explist is None:
                continue
            for exp in explist:
                if exp.fullmatch(filepath):
                    return cmd
        return None
