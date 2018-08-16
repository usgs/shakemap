.. _sec-output-formats-4:

****************************
Output Data Formats
****************************

About the HDF5 File Format
==========================

HDF5 files are an efficient and flexible way to store large sets of
structured data.  HDF stands for Hierarchical Data Format and is an
open-source project . The data model consists of two types of objects:
groups and datasets. Groups may contain datasets or other groups
(along with associated metadata).  Datasets may contain
multidimensional arrays of data (along with accompanying
metadata). Groups may be thought of as directories or folders, and
datasets as files. The metadata associated with each object allows
HDF5 files to be self describing. Data is stored in binary, making the
files smaller and more efficient to read and write than text-based
formats (such as XML). For almost all uses a user will access the file
through an API (bindings are available for many languages), making the
details of file structure and layout not relevant.

For more information on HDF5 see the 
`HDF5 support page <https://support.hdfgroup.org/HDF5/>`_.
The Python interface, h5py, is documented on the 
`h5py web site <http://www.h5py.org/>`_.

For ShakeMap, the two main data files (*shake_data.hdf* and 
*shake_result.hdf*) are accessed through the methods of the
shakelib classes **ShakeMapInputContainer** and **ShakeMapOutputContainer**.
For documentation see the `shakelib container module 
<https://usgs.github.io/shakelib/shakelib.utils.container.html>`_.

shake_data.hdf
================

*shake_data.hdf* is the input to the **model** module. It is built by 
bringing 
together the various ShakeMap configuration files (found in both the 
current profile's *install/config* path, and in the event's *current*
directory. It also aggregates information from the event's rupture file
(if one exists), *event.xml*, and the data files -- all found in the 
event's *current* directory. **shake_data.hdf** is usually created 
with the **assemble** module (and possibly modified by **augment**). 

It is generally not necessary for operators to access this file other
than through the interfaces of **assemble** and **augment**. For
developers, the file
is accessed through the shakelib ShakeMapInputContainer interface:
https://usgs.github.io/shakemap/shakelib/shakelib.utils.containers.html

shake_result.hdf
================

The primary output product of the **model** module is the HDF5 file 
*shake_result.hdf*. There are some differences depending on whether 
**model** produces output grids or lists of points.  HDF5 files
consist of *Groups* or *DataSets*.  Groups can contain other groups
or datasets, and can contain a dictionary-like set of *attributes*.
Datasets consist of arrays holding data values, and
metadata in the form of attributes.  A very good introduction to
the HDF5 data format can be found
here: https://support.hdfgroup.org/HDF5/Tutor/HDF5Intro.pdf

*shake_result.hdf* consists of a number of groups and datasets. Our
implementation of HDF5 uses groups to contain Python dictionaries,
strings, and numpy arrays.  Dictionaries are stored as JSON strings.
For example, a Python *config* dictionary consisting of the following information:

.. code-block:: python

    d = {
        'name': 'config',
        'gmpe_set': {
            'Active_Crustal': 'Campbell2003'
        },
        'depth': 34.0,
    }

would be stored as a JSON string in a dataset called ``config``,
under a group called ``dictionaries``. In languages that support
JSON, these strings can be easily converted into a data
structure. Developers can access these properties using the
ShakeMapInputContainer interface:
https://usgs.github.io/shakemap/shakelib/shakelib.utils.containers.html

*shake_result.hdf* contains the following metadata elements:

+-----------------------+----------------------------+-------------+----------------------------------------------+
| Name                  | Location                   | Python Type | Contents                                     |
+=======================+============================+=============+==============================================+
| config                | /dictionaries/config       | JSON string | ShakeMap configuration                       |
+-----------------------+----------------------------+-------------+----------------------------------------------+
| info.json             | /dictionaries/info.json    | JSON string | ShakeMap metadata                            |
+-----------------------+----------------------------+-------------+----------------------------------------------+
| stations_dict         | /dictionaries/stations_dict| JSON string | Dictionary representation of observed data   |
+-----------------------+----------------------------+-------------+----------------------------------------------+
| rupture               | /dictionaries/rupture      | JSON string | Dictionary representation of fault rupture   |
+-----------------------+----------------------------+-------------+----------------------------------------------+

It also will contain a number of grids, which, when read with the HDFContainer getGrid() method, return
a Grid2D object which is a Python representation of a North-up 2D array of data, whose upper-left corner
coordinate and cell dimensions are known.  The definition of this object can be found here:
https://github.com/usgs/MapIO/blob/master/mapio/grid2d.py


Sampling of grids contained in the HDF:

+-----------------------+----------------------------+-------------+----------------------------------------------+
| Name                  | Location                   | Python Type | Contents                                     |
+=======================+============================+=============+==============================================+
| vs30                  | /grids/vs30                | Grid2D      | Vs30 values at output grid/points            | 
+-----------------------+----------------------------+-------------+----------------------------------------------+
| distance_rhypo        | /grids/distance_rhypo      | Grid2D      | Hypocentral distance                         |
+-----------------------+----------------------------+-------------+----------------------------------------------+
| distance_rjb          | /grids/distance_rjb        | Grid2D      | RJB distance                                 |
+-----------------------+----------------------------+-------------+----------------------------------------------+


Each IMT dataset (MMI, PGA, etc.) is stored as a group containing two 
datasets, the mean values for each cell and the standard deviations.  
MMI data for the component 'Larger' will be stored under a group called 
``imt_MMI_GREATER_OF_TWO_HORIZONTAL``. The mean array will be stored as
``mean``, and the standard deviation array will be stored as
``std``.  All IMT grid datasets will be accompanied by a dictionary of
attributes:

+-----------+------------------------------------------------------+
| Attr name | Contents                                             |
+===========+======================================================+
| units     | Physical units of IMT or standard deviation.         |
+-----------+------------------------------------------------------+
| digits    | Number of significant digits to use for the values.  |
+-----------+------------------------------------------------------+
| xmin      | The eastern boundary of the grid (degrees longitude) |
+-----------+------------------------------------------------------+
| xmax      | The eastern boundary of the grid (degrees longitude) |
+-----------+------------------------------------------------------+
| ymin      | The southern boundary of the grid (degrees latitude) |
+-----------+------------------------------------------------------+
| ymax      | The northern boundary of the grid (degrees latitude) |
+-----------+------------------------------------------------------+
| nx        | The number of grid points in the x dimension         |
+-----------+------------------------------------------------------+
| ny        | The number of grid points in the y dimension         |
+-----------+------------------------------------------------------+
| dx        | The grid interval in the x dimension                 |
+-----------+------------------------------------------------------+
| dy        | The grid interval in the y dimension                 |
+-----------+------------------------------------------------------+

Sampling of IMTs in the HDF file:

+-----------------------+----------------------------------------------+----------------------------------------------+
| Name                  | Location                                     | Python Type | Contents                       |
+=======================+==============================================+=============+================================+
| MMI Mean              | /imts/MMI_GREATER_OF_TWO_HORIZONTAL/mean     | Grid2D      | MMI Mean Values                | 
+-----------------------+----------------------------+-----------------+----------------------------------------------+
| MMI Std               | /imts/MMI_GREATER_OF_TWO_HORIZONTAL/std      | Grid2D      | MMI Std                        | 
+-----------------------+----------------------------+-----------------+----------------------------------------------+
| Sa(0.3) Mean          | /imts/SA(0.3)_GREATER_OF_TWO_HORIZONTAL/mean | Grid2D      | SA(0.3) Mean Values            | 
+-----------------------+----------------------------+-----------------+----------------------------------------------+
| Sa(0.3) Std           | /imts/SA(0.3)_GREATER_OF_TWO_HORIZONTAL/std  | Grid2D      | SA(0.3) Std                    | 
+-----------------------+----------------------------+-------------+--------------------------------------------------+

For datasets that are lists of points, the storage of IMTS is the same
as for grids, except that the data are stored as one-dimensional arrays.
Each IMT group wll also contains datasets ``lons``, ``lats``, 
and ``ids``, which provide the coordinates of the points in longitude
and latitude, and their IDs, respectively. For sets of points the metadata
attributes are:

+--------------+------------------------------------------------------+
| Attr name    | Contents                                             |
+==============+======================================================+
| units        | Physical units of the IMT                            |
+--------------+------------------------------------------------------+
| digits       | Number of significant digits to use for the values   |
+--------------+------------------------------------------------------+

All *shake_result.hdf* files will have a group ``file_data_type`` 
which will have a single attribute ``data_type`` that will be one of
'points' or 'grid'. This way the user can distinguish between the two
types of storage.

For grid files, there will also be sets of regression curves stored
as one-dimensional arrays. These
will be labeled like ``array_regression_<IMT>_<site>_<type>`` Where
``<IMT>`` will be one of the output IMTs (e.g., ``SA(3.0)``), 
``<site>`` will be one of ``rock`` or ``soil`` (for which ``rock``
means a Vs30 of 760 m/s, and ``soil`` means a Vs30 of 180 m/s), and
``<type>`` is one of ``mean`` (for the mean values) or ``sd`` (for
the standard deviations). All units are in natural log space (except
for MMI). There will also be an array called 
``array_regression_distances`` which will contain the distances
(in km) corresponding to the points in the data arrays.

Regardless of whether the file stores grids or arrays of points, it will
also contain datasets of various distance parameters. These will be 
named ``distance_*`` where the wildcard will be replaced with  one of 
the typical 
source distance metrics (e.g., 'rrup' for rupture distance, 'rjb' for
Joyner-Boore distance, 'rhypo' for hypocentral distance, 'rx', 'ry0',
etc.) 
Similarly, the Vs30 data are found in a dataset named ``vs30``.
The metadata for distances and Vs30 consists of 'units' and 'digits'.

Dictionary datasets are stored as JSON strings.

There will typically be multiple *IMT* (Intensity Measure Type) datasets
(each containing the mean and standard deviation of the IMT). For instance
'PGA', 'PGV', 'MMI', and various 'SA(#num)' 
[where #num is the period as a floating point number; e.g., 
'SA(1.0)']. 

Python developers will likely want to access *shake_result.hdf* through
the shakelib OutputContainer class:
https://usgs.github.io/shakelib/shakelib.utils.container.html
Also see, for example, the *contour* module [:meth:`shakemap.coremods.contour`]
for some basic access patterns.

Matlab developers can use the function *read_shake_data.m*, which is included in
the repository for ShakeMap here: https://github.com/usgs/shakemap/blob/master/contrib/read_shakemap_data.m


