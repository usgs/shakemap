.. _sec-formats-4:

****************************
Data Formats
****************************

About the HDF5 File Format
==========================

HDF5 files are an efficient and flexible way to store large sets of 
structured data.  HDF stands for Hierachical Data Format and is an 
open-source project . The data 
model consists of two types of objects: groups and datasets. Groups
may contain datasets or other groups (along with associated metadata).
Datasets may contain multidimensional arrays of data (along with
accompanying metadata). Groups may be thought of as directories or
folders, and datasets as files. The metadata associated with each object
allows HDF5 files to be self describing. Data is stored in binary,
making the files smaller and more efficient to read and write than 
text-based formats (such as XML). For almost all uses
a user will access the file through an API (bindings are available
for many languages), making the details of file structure and layout 
not relevant.

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
is accessed through the shakelib `ShakeMapInputContainer interface 
<https://usgs.github.io/shakelib/shakelib.utils.container.html>`_.

shake_result.hdf
================

The primary output product of the **model** module is the HDF5 file 
*shake_result.hdf*. There are some differences depending on whether 
**model** produces output grids or lists of points.  HDF5 files
consist of *Groups* or *DataSets*.  Groups can contain other groups
or datasets, and can contain a dictionary-like set of *attributes*.
Datasets consist of arrays holding data values, and
metadata in the form of attributes.  A very good introduction to
the HDF5 data format can be found here:

https://support.hdfgroup.org/HDF5/Tutor/HDF5Intro.pdf

*shake_result.hdf* consists of a number of groups and datasets. Our
implementation of HDF5 uses groups to contain Python dictionaries,
strings, and numpy arrays.  Dictionaries are stored as recursive groups.
For example, a Python *config* dictionary consisting of the following information::

  d = {
       'name': 'config',
       'gmpe_set': {
           'Active_Crustal': 'Campbell2003'
        },
       'depth': 34.0,
  }

would be stored in a group called ``__dictionary_config__``, with attributes
'name' and 'depth'.  There will be a sub-group called 'gmpe_set', with an attribute
called 'Active_Crustal'.  Groups like this can contain any number (within reason)
of dictionaries, which can consist of a number of Python data types (strings, numbers,
numpy arrays, etc.) *shake_result.hdf* contains the following elements:  

+-----------------------+---------+-------------+-----------------------------------------+
| Name                  | HDF Type| Python Type | Contents                                |
+=======================+=========+=============+=========================================+
| config                | group   | dictionary  | ShakeMap configuration                  |
+-----------------------+---------+-------------+-----------------------------------------+
| info.json             | group   | string      | ShakeMap metadata in JSON format        |
+-----------------------+---------+-------------+-----------------------------------------+
| rupture.json          | group   | string      | GeoJSON representation of the           |
|                       |         |             | finite fault                            |
+-----------------------+---------+-------------+-----------------------------------------+
| stationlist.json      | dataset | string      | GeoJSON stationlist                     |
+-----------------------+---------+-------------+----------------------+------------------+
| vs30                  | dataset | ndarray     | Vs30 values at output grid/points       | 
+-----------------------+---------+-------------+----------------------+------------------+
| IMT (multiple)        | group   | dictionary  | Interpolated data for IMT               |
+-----------------------+---------+-------------+----------------------+------------------+

Each IMT dataset (MMI,PGA,etc.) is stored as a group containing two datasets, the mean values
for each cell and the standard deviations.  MMI data for the component 'Larger' will be stored
under a group called ``__imt_MMI_Larger__``. The mean array will be stored as
``__mean_MMI_Larger__``, and the standard deviation array will be stored as
``__std_MMI_Larger__``.  All IMT grid datasets will be accompanied by a dictionary of
attributes:


For grids, the attributes are:

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

For lists of points the attributes are:

+--------------+------------------------------------------------------+
| Attr name    | Contents                                             |
+==============+======================================================+
| type         | "points"                                             |
+--------------+------------------------------------------------------+
| lons         | An array of longitudes corresponding to the points   |
+--------------+------------------------------------------------------+
| lats         | An array of latitudes corresponding to the points    |
+--------------+------------------------------------------------------+
| facility_ids | array of identifiers for the points                  |
+--------------+------------------------------------------------------+

JSON datasets are stored as continuous strings.

There will typically be multiple *IMT* (Intensity Measure Type) and 
*IMT_sd* (standard deviation of the corresponding IMT) datasets. For instance
*PGA*, *PGA_sd*, *PGV*, *PGV_sd*, *MMI*, *MMI_sd*, and various *SA(#num)* and
*SA(#num)_sd* [where #num is the period as a floating point number; e.g., 
*SA(1.0)*]. 

Python developers will likely want to access *shake_result.hdf* through
the `shakelib OutputContainer class 
<https://usgs.github.io/shakelib/shakelib.utils.container.html>`_.
Also see, for example, the *contour* module [:meth:`shakemap.coremods.contour`]
for some basic access patterns.

Matlab developers can use the function *read_shake_data.m*, which is included in
the repository for ShakeMap
(https://github.com/usgs/shakemap/blob/master/read_shakemap_data.m).
