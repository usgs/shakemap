#!/usr/bin/env python

# stdlib imports
import json
import logging

# third party imports
from mapio.grid2d import Grid2D
from mapio.geodict import GeoDict
from mapio.gridcontainer import GridHDFContainer, _split_dset_attrs
from impactutils.io.container import _get_type_list

# local imports
from shakelib.rupture.factory import (rupture_from_dict,
                                      get_rupture,
                                      rupture_from_dict_and_origin)
from shakelib.rupture.origin import Origin
from shakelib.station import StationList

class ShakeMapContainer(GridHDFContainer):
    """
    Parent class for InputShakeMapContainer and OutputShakeMapContainer.
    """
    def setConfig(self, config):
        """
        Add the config as a dictionary to the HDF file.

        Args:
            config (dict--like): Dict--like object with configuration
                information.
        """
        if 'config' in self.getDictionaries():
            self.dropDictionary('config')
        self.setDictionary('config', config)

    def getConfig(self):
        """
        Retrieve configuration dictionary from container.

        Returns:
            dict: Configuration dictionary.
        Raises:
            AttributeError: If config dictionary has not been set in
                the container.
        """
        if 'config' not in self.getDictionaries():
            raise AttributeError('Configuration not set in container.')
        config = self.getDictionary('config')
        return config

    def setRupture(self,rupture):
        """
        Store Rupture object in container.

        Args:
            rupture (dict or Rupture): Rupture object (Point,Quad, or Edge)
                or dictionary representation of same.
        Raises:
            TypeError: If input object or dictionary does not
                represent a Rupture object.
        """
        if 'rupture' in self.getStrings():
            self.dropString('rupture')
        if isinstance(rupture,dict):
            try:
                _ = rupture_from_dict(rupture)
            except Exception:
                fmt = 'Input dict does not represent a rupture object.'
                raise TypeError(fmt)
            json_str = json.dumps(rupture)
        else:
            json_str = json.dumps(rupture._geojson)
        self.setString('rupture',json_str)

    def getRuptureObject(self):
        """
        Retrieve Rupture object from container.

        Returns:
            Rupture: Instance of (one of) a Point/Quad/EdgeRupture class.
        Raises:
            AttributeError: If rupture object has not been set in
                the container.
        """
        rupture_dict = self.getRuptureDict()
        rupture = rupture_from_dict(rupture_dict)
        return rupture

    def getRuptureDict(self):
        """
        Retrieve Rupture dictionary from container.

        Returns:
            dict: Dictionary representatin of (one of) a
                Point/Quad/EdgeRupture class.
        Raises:
            AttributeError: If rupture object has not been set in
                the container.
        """
        if 'rupture' not in self.getStrings():
            raise AttributeError('Rupture object not set in container.')
        rupture_dict = json.loads(self.getString('rupture'))
        return rupture_dict

    def setStationList(self,stationlist):
        """
        Store StationList object in container.

        Args:
            stationlist (StationList): StationList object.
        Raises:
            TypeError: If input object or dictionary is not a StationList object.
        """
        if 'stations' in self.getStrings():
            self.dropString('stations')
        if not isinstance(stationlist,StationList):
            fmt = 'Input object is not a StationList.'
            raise TypeError(fmt)
        sql_string = stationlist.dumpToSQL()
        self.setString('stations',sql_string)

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
        sql_string = self.getString('stations')
        stationlist = StationList.loadFromSQL(sql_string)
        return stationlist

    def setStationDict(self, stationdict):
        """
        Store (JSON-like) station dictiionary in container.

        Args:
            stationdict (dict-like): Station dict object.
        Raises:
            TypeError: If input object is not a dictionary.
        """
        if not isinstance(stationdict, dict):
            fmt = 'Input object is not a dictionary.'
            raise TypeError(fmt)
        if 'stations_dict' in self.getStrings():
            self.dropStrings('stations_dict')
        stationstring = json.dumps(stationdict)
        self.setString('stations_dict', stationstring)

    def getStationDict(self):
        """
        Retrieve (JSON-like) station dictionary from container.

        Returns:
            dict-like: Station dictionary.
        Raises:
            AttributeError: If station dictionary has not been set in
                the container.
        """
        if 'stations_dict' not in self.getStrings():
            raise AttributeError('Station dictionary not set in container.')
        stations = json.loads(self.getString('stations_dict'))
        return stations

    def setVersionHistory(self, history_dict):
        """
        Store a dictionary containing version history in the container.

        Args:
            history_dict (dict): Dictionary containing version history. ??
        """
        if 'version_history' in self.getDictionaries():
            self.dropDictionary('version_history')
        self.setDictionary('version_history', history_dict)

    def getVersionHistory(self):
        """
        Return the dictionary containing version history.

        Returns:
          dict: Dictionary containing version history, or None.
        Raises:
            AttributeError: If version history has not been set in
                the container.
        """

        if 'version_history' not in self.getDictionaries():
            raise AttributeError('Version history not set in container.')

        version_dict = self.getDictionary('version_history')
        return version_dict


class ShakeMapInputContainer(ShakeMapContainer):
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
    def createFromInput(cls, filename, config, eventfile, version_history, rupturefile=None,
                        datafiles=None):
        """
        Instantiate an InputContainer from ShakeMap input data.

        Args:
            filename (str): Path to HDF5 file that will be created to encapsulate all
                input data.
            config (dict): Dictionary containing all configuration information
                necessary for ShakeMap ground motion and other calculations.
            eventfile (str): Path to ShakeMap event.xml file.
            rupturefile (str): Path to ShakeMap rupture text or JSON file.
            datafiles (list): List of ShakeMap data (DYFI, strong motion) files.
            version_history (dict): Dictionary containing version history.

        Returns:
            InputContainer: Instance of InputContainer.
        """
        container = cls.create(filename)
        container.setConfig(config)

        # create an origin from the event file
        origin = Origin.fromFile(eventfile)

        # create a rupture object from the origin and the rupture file
        # (if present).
        rupture = get_rupture(origin,file=rupturefile)

        # store the rupture object in the container
        container.setRupture(rupture)

        if datafiles is not None:
            container.setStationData(datafiles)

        container.setVersionHistory(version_history)
        return container

    def setStationData(self, datafiles):
        """
        Insert observed ground motion data into the container.

        Args:
          datafiles (str): Path to an XML-formatted file containing ground
              motion observations, (macroseismic or instrumented).

        """
        station = StationList.loadFromXML(datafiles)
        self.setStationList(station)

    def addStationData(self, datafiles):
        """
        Add observed ground motion data into the container.

        Args:
            datafiles (str): Path to an XML-formatted file containing ground
                motion observations, (macroseismic or instrumented).

        """
        station = self.getStationList()
        station.addData(datafiles)
        self.setStationList(station)

    def updateRupture(self,eventxml = None, rupturefile = None):
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
                rupture = get_rupture(origin,file=rupturefile)
            else:
                rupture_dict = rupture._geojson
                rupture = rupture_from_dict_and_origin(rupture_dict,origin)
        else: #no event.xml file, but we do have a rupture file
            rupture = get_rupture(origin,file=rupturefile)

        self.setRupture(rupture)


class ShakeMapOutputContainer(ShakeMapContainer):
    """
    HDF container for Shakemap output data.

    This class provides methods for getting and setting IMT data.
    The philosophy here is that an IMT consists of both the mean results and
    the standard deviations of those results, thus getIMT() returns a
    dictionary with both, plus metadata for each data layer.


    """

    def getDataType(self):
        """
        Returns the format of the IMT, Vs30, and distance data stored in
        this file: either 'points' or 'grid'. None is returned if no data
        have been set.

        Returns:
            str or None: Either 'grid' or 'points' or None.
        """

        group_name = '__file_data_type__'
        if group_name in self._hdfobj:
            data_type_group = self._hdfobj[group_name]
            return data_type_group.attrs['data_type']
        return None

    def setDataType(self, datatype):
        """
        Sets the type of the IMT, Vs30, and distance data stored in this
        file. This function should not be called by the user -- the value
        will be set by the first call to setIMTGrids() or setIMTArray().

        Args:
            (str): Either 'points' or 'grid'.

        Returns:
            Nothing.
        """

        if datatype != 'points' and datatype != 'grid':
            raise TypeError('Trying to set unknown data type: %s' %
                            (datatype))
        group_name = '__file_data_type__'
        if group_name in self._hdfobj:
            data_type_group = self._hdfobj[group_name]
            current_data_type = data_type_group.attrs['data_type']
            if current_data_type != datatype:
                raise TypeError(
                        'Trying to set data type to %s; file already type %s'
                        % (datatype, current_data_type))
            #
            # Data type is already set; don't have to do anything
            #
            return
        data_type_group = self._hdfobj.create_group(group_name)
        data_type_group.attrs['data_type'] = datatype
        return

    def setIMTGrids(self, imt_name, imt_mean, mean_metadata,
                    imt_std, std_metadata, component,
                    compression=True):
        """
        Store IMT mean and standard deviation objects as datasets.

        Args:
            imt_name (str): Name of the IMT (MMI, PGA, etc.) to be stored.
            imt_mean (Grid2D): Grid2D object of IMT mean values to be stored.
            mean_metadata (dict): Dictionary containing metadata for mean IMT
                grid.
            imt_std (Grid2D): Grid2D object of IMT standard deviation values
                to be stored.
            std_metadata (dict): Dictionary containing metadata for mean IMT
                grid.
            component (str): Component type, i.e. 'Larger','rotd50',etc.
            compression (bool): Boolean indicating whether dataset should be
                compressed using the gzip algorithm.

        Returns:
            HDF Group containing IMT grids and metadata.
        """

        if self.getDataType() == 'points':
            raise TypeError('Setting grid data in a file containing points')
        self.setDataType('grid')

        # set up the name of the group holding all the information for the IMT
        group_name = '__imt_%s_%s__' % (imt_name, component)
        if group_name in self._hdfobj:
            raise Exception('An IMT group called %s already exists.' % imt_name)
        imt_group = self._hdfobj.create_group(group_name)

        # create the data set containing the mean IMT data and metadata
        mean_data = imt_mean.getData()
        mean_set = imt_group.create_dataset('mean', data=mean_data,
                                            compression=compression)
        if mean_metadata is not None:
            for key, value in mean_metadata.items():
                mean_set.attrs[key] = value
        for key, value in imt_mean.getGeoDict().asDict().items():
            mean_set.attrs[key] = value

        # create the data set containing the std IMT data and metadata
        std_data = imt_std.getData()
        std_set = imt_group.create_dataset('std', data=std_data,
                                           compression=compression)
        if std_metadata is not None:
            for key, value in std_metadata.items():
                std_set.attrs[key] = value
        for key, value in imt_std.getGeoDict().asDict().items():
            std_set.attrs[key] = value

        return imt_group

    def getIMTGrids(self, imt_name, component):
        """
        Retrieve a Grid2D object and any associated metadata from the container.

        Args:
            imt_name (str):
                The name of the IMT stored in the container.

        Returns:
            dict: Dictionary containing 4 items:
                   - mean Grid2D object for IMT mean values.
                   - mean_metadata Dictionary containing any metadata
                     describing mean layer.
                   - std Grid2D object for IMT standard deviation values.
                   - std_metadata Dictionary containing any metadata describing
                     standard deviation layer.
        """

        if self.getDataType() != 'grid':
            raise TypeError('Requesting grid data from file containing points')

        group_name = '__imt_%s_%s__' % (imt_name, component)
        if group_name not in self._hdfobj:
            raise LookupError('No group called %s in HDF file %s'
                              % (imt_name, self.getFileName()))
        imt_group = self._hdfobj[group_name]

        # get the mean data and metadata
        mean_dset = imt_group['mean']
        mean_data = mean_dset[()]

        array_metadata, mean_metadata = _split_dset_attrs(mean_dset)
        mean_geodict = GeoDict(array_metadata)
        mean_grid = Grid2D(mean_data, mean_geodict)

        # get the std data and metadata
        std_dset = imt_group['std']
        std_data = std_dset[()]

        array_metadata, std_metadata = _split_dset_attrs(std_dset)
        std_geodict = GeoDict(array_metadata)
        std_grid = Grid2D(std_data, std_geodict)

        # create an output dictionary
        imt_dict = {
            'mean': mean_grid,
            'mean_metadata': mean_metadata,
            'std': std_grid,
            'std_metadata': std_metadata
        }
        return imt_dict

    def setIMTArrays(self, imt_name, lons, lats, ids,
                     imt_mean, mean_metadata,
                     imt_std, std_metadata,
                     component, compression=True):
        """
        Store IMT mean and standard deviation objects as datasets.

        Args:
            imt_name (str): Name of the IMT (MMI, PGA, etc.) to be stored.
            lons (Numpy array): Array of longitudes of the IMT data.
            lats (Numpy array): Array of latitudes of the IMT data.
            ids (Numpy array): Array of ID strings corresponding to the
                locations given by lons and lats.
            imt_mean (Numpy array): Array of IMT mean values to be stored.
            mean_metadata (dict): Dictionary containing metadata for mean IMT
                grid.
            imt_std (Numpy array): Array of IMT standard deviation values
                to be stored.
            std_metadata (dict): Dictionary containing metadata for mean IMT
                grid.
            component (str): Component type, i.e. 'Larger','rotd50',etc.
            compression (bool): Boolean indicating whether dataset should be
                compressed using the gzip algorithm.

        Returns:
            HDF Group containing IMT arrays and metadata.
        """

        if self.getDataType() == 'grid':
            raise TypeError('Setting point data in a file containing grids')
        self.setDataType('points')

        #
        # Check that all of the arrays are the same
        # size
        #
        if lons.shape != lats.shape or \
           lons.shape != ids.shape or \
           lons.shape != imt_mean.shape or \
           lons.shape != imt_std.shape:
            raise ValueError('All input arrays must be the same shape')

        # set up the name of the group holding all the information for the IMT
        group_name = '__imt_%s_%s__' % (imt_name, component)
        if group_name in self._hdfobj:
            raise ValueError('An IMT group called %s already exists.' %
                            imt_name)
        imt_group = self._hdfobj.create_group(group_name)

        # create data sets containing the longitudes, latitudes, and ids
        imt_group.create_dataset('lons', data=lons, compression=compression)
        imt_group.create_dataset('lats', data=lats, compression=compression)
        imt_group.create_dataset('ids', data=ids, compression=compression)

        # create the data set containing the mean IMT data and metadata
        dataset = imt_group.create_dataset('mean', data=imt_mean,
                                           compression=compression)
        if mean_metadata is not None:
            for key, value in mean_metadata.items():
                dataset.attrs[key] = value

        # create the data set containing the std IMT data and metadata
        dataset = imt_group.create_dataset('std', data=imt_std,
                                           compression=compression)
        if std_metadata is not None:
            for key, value in std_metadata.items():
                dataset.attrs[key] = value

        return imt_group

    def getIMTArrays(self, imt_name, component):
        """
        Retrieve the arrays and any associated metadata from the container.

        Args:
            imt_name (str):
                The name of the IMT stored in the container.

        Returns:
            dict: Dictionary containing 7 items:
                   - 'lons': array of longitude coordinates
                   - 'lats': array of latitude coordinates
                   - 'ids': array of IDs corresponding to the coordinates
                   - 'mean': array of IMT mean values.
                   - 'mean_metadata': Dictionary containing any metadata
                     describing mean layer.
                   - 'std': array of IMT standard deviation values.
                   - 'std_metadata': Dictionary containing any metadata
                     describing standard deviation layer.
        """

        if self.getDataType() != 'points':
            raise TypeError('Requesting point data from file containing grids')

        group_name = '__imt_%s_%s__' % (imt_name, component)
        if group_name not in self._hdfobj:
            raise LookupError('No group called %s in HDF file %s'
                              % (imt_name, self.getFileName()))
        imt_group = self._hdfobj[group_name]

        # get the coordinates and ids
        dset = imt_group['lons']
        lons = dset[()]
        dset = imt_group['lats']
        lats = dset[()]
        dset = imt_group['ids']
        ids = dset[()]

        # get the mean data and metadata
        mean_dset = imt_group['mean']
        mean_data = mean_dset[()]

        mean_metadata = {}
        for key in mean_dset.attrs.keys():
            mean_metadata[key] = mean_dset.attrs[key]

        # get the std data and metadata
        std_dset = imt_group['std']
        std_data = std_dset[()]

        std_metadata = {}
        for key in std_dset.attrs.keys():
            std_metadata[key] = std_dset.attrs[key]

        # create an output dictionary
        imt_dict = {
            'lons': lons,
            'lats': lats,
            'ids': ids,
            'mean': mean_data,
            'mean_metadata': mean_metadata,
            'std': std_data,
            'std_metadata': std_metadata
        }
        return imt_dict

    def getIMTs(self, component):
        """
        Return list of names of IMTs matching input component type.

        Args:
            component (str): Name of component ('Larger','rotd50',etc.)

        Returns:
            list: List of names of IMTs matching component stored in container.
        """
        imt_groups = _get_type_list(self._hdfobj, 'imt')
        comp_groups = []
        for imt_group in imt_groups:
            if imt_group.find(component) > -1:
                comp_groups.append(imt_group.replace('_'+component, ''))
        return comp_groups

    def getComponents(self, imt_name):
        """
        Return list of components for given IMT.

        Args:
          imt_name (str): Name of IMT ('MMI','PGA',etc.)

        Returns:
          list: List of names of components for given IMT.
        """
        components = _get_type_list(self._hdfobj, 'imt_'+imt_name)
        return components

    def dropIMT(self, imt_name):
        """
        Delete IMT datasets from container.

        Args:
          name (str):
                The name of the IMT to be deleted.

        """
        group_name = '__imt_%s_%s__' % (imt_name, component)
        if group_name not in self._hdfobj:
            raise LookupError('No group called %s in HDF file %s'
                              % (imt_name, self.getFileName()))
        del self._hdfobj[group_name]
