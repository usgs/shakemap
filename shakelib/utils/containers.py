#!/usr/bin/env python

# stdlib imports

# third party imports
from impactutils.io.smcontainers import ShakeMapContainerBase

# local imports
from shakelib.rupture.factory import (rupture_from_dict,
                                      get_rupture,
                                      rupture_from_dict_and_origin)
from shakelib.rupture.origin import Origin
from shakelib.station import StationList


class ShakeMapInputContainer(ShakeMapContainerBase):
    """
    HDF container for Shakemap input data.

    This class provides methods for getting and setting information on:
     - configuration
     - event (lat,lon,depth,etc.)
     - rupture
     - station data (can also be updated)
     - version history

    """
    @classmethod
    def createFromInput(cls, filename, config, eventfile, version_history,
                        rupturefile=None, sourcefile=None, momentfile=None,
                        datafiles=None):
        """
        Instantiate an InputContainer from ShakeMap input data.

        Args:
            filename (str): Path to HDF5 file that will be created to
                encapsulate all input data.
            config (dict): Dictionary containing all configuration information
                necessary for ShakeMap ground motion and other calculations.
            eventfile (str): Path to ShakeMap event.xml file.
            rupturefile (str): Path to ShakeMap rupture text or JSON file.
            sourcefile (str): Path to ShakeMap source.txt file.
            momentfile (str): Path to ShakeMap moment.xml file.
            datafiles (list): List of ShakeMap data (DYFI, strong motion)
                files.
            version_history (dict): Dictionary containing version history.

        Returns:
            InputContainer: Instance of InputContainer.
        """
        container = cls.create(filename)
        container.setConfig(config)

        # create an origin from the event file
        origin = Origin.fromFile(eventfile, sourcefile=sourcefile,
                                 momentfile=momentfile)

        # create a rupture object from the origin and the rupture file
        # (if present).
        rupture = get_rupture(origin, file=rupturefile)

        # store the rupture object in the container
        container.setRupture(rupture)

        if datafiles is not None and len(datafiles) > 0:
            container.setStationData(datafiles)

        container.setVersionHistory(version_history)
        return container

    def setRupture(self, rupture):
        """
        Store Rupture object in container.

        Args:
            rupture (dict or Rupture): Rupture object (Point,Quad, or Edge)
                or dictionary representation of same.
        Raises:
            TypeError: If input object or dictionary does not
                represent a Rupture object.
        """
        if 'rupture' in self.getDictionaries():
            self.dropDictionary([], 'rupture')
        if isinstance(rupture, dict):
            try:
                rupture_from_dict(rupture)
            except Exception:
                fmt = 'Input dict does not represent a rupture object.'
                raise TypeError(fmt)
            rupdict = rupture
        else:
            rupdict = rupture._geojson
        self.setRuptureDict(rupdict)

    def getRuptureObject(self):
        """
        Retrieve Rupture object from container.

        Returns:
            Rupture: Instance of (one of) a Point/Quad/EdgeRupture class.
        Raises:
            AttributeError: If rupture object has not been set in
                the container.
        """
        if 'rupture' not in self.getDictionaries():
            raise AttributeError('Rupture dictionary not set in container.')
        rupture_dict = self.getRuptureDict()
        rupture = rupture_from_dict(rupture_dict)
        return rupture

    def setStationList(self, stationlist):
        """
        Store StationList object in container.

        Args:
            stationlist (StationList): StationList object.
        Raises:
            TypeError: If input object or dictionary is not a StationList
                object.
        """
        if 'stations' in self.getStrings():
            self.dropString([], 'stations')
        if not isinstance(stationlist, StationList):
            fmt = 'Input object is not a StationList.'
            raise TypeError(fmt)
        sql_string = stationlist.dumpToSQL()
        self.setString([], 'stations', sql_string)

    def getStationList(self):
        """
        Retrieve StationList object from container.

        Returns:
            StationList: StationList object.
        Raises:
            AttributeError: If stationlist object has not been set in
                the container.
        """
        if 'stations' not in self.getStrings():
            raise AttributeError('StationList object not set in container.')
        sql_string = self.getString([], 'stations')
        stationlist = StationList.loadFromSQL(sql_string)
        return stationlist

    def setStationData(self, datafiles):
        """
        Insert observed ground motion data into the container.

        Args:
          datafiles (str): Path to XML- or JSON-formatted files containing
              ground motion observations, (macroseismic or instrumented).

        """
        station = StationList.loadFromFiles(datafiles)
        self.setStationList(station)

    def addStationData(self, datafiles):
        """
        Add observed ground motion data into the container.

        Args:
            datafiles (sequence): Sequence of paths to XML- and/or
                JSON-formatted files containing ground motion observations,
                (macroseismic or instrumented).
        """
        if len(datafiles) == 0:
            return
        try:
            station = self.getStationList()
            station.addData(datafiles)
        except AttributeError:
            station = StationList.loadFromFiles(datafiles)
        self.setStationList(station)

    def updateRupture(self, eventxml=None, rupturefile=None):
        """Update rupture/origin information in container.

        Args:
            eventxml (str): Path to event.xml file.
            rupturefile (str): Path to rupture file (JSON or txt format).
        """
        if eventxml is None and rupturefile is None:
            return

        # the container is guaranteed to have at least a Point rupture
        # and the origin.
        rupture = self.getRuptureObject()
        origin = rupture.getOrigin()

        if eventxml is not None:
            origin = Origin.fromFile(eventxml)
            if rupturefile is not None:
                rupture = get_rupture(origin, file=rupturefile)
            else:
                rupture_dict = rupture._geojson
                rupture = rupture_from_dict_and_origin(rupture_dict, origin)
        else:  # no event.xml file, but we do have a rupture file
            rupture = get_rupture(origin, file=rupturefile)

        self.setRupture(rupture)
