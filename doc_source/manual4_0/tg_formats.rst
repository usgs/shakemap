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
shakelib classes **InputContainer** and **OutputContainer**.
For documentation see the `shakelib container module 
<https://usgs.github.io/shakelib/shakelib.utils.container.html>`_.

shake_data.hdf
================

**shake_data.hdf** is the input to *sm_model*. It is built by bringing 
together the various ShakeMap configuration files (found in both the 
current profile's ``install/config`` path, and in the event's ``current``
directory. It also aggregates information from the event's rupture file
(if one exists), ``event.xml``, and the data files -- all found in the 
event's ``current`` directory. **shake_data.hdf** is usually created 
with *sm_assemble* (and possibly modified by *sm_augment*). 

It is generally not necessary for operators to access this file other
than through the interfaces of *sm_assemble* and *sm_augment*. For
developers, the file
is accessed through the shakelib `InputContainer interface 
<https://usgs.github.io/shakelib/shakelib.utils.container.html>`_.

shake_result.hdf
================

The primary output product of **sm_model** is the HDF5 file 
*shake_result.hdf*. There are some differences depending on whether 
**sm_model** produces output grids or lists of points. *shake_result.hdf*
is essentially "flat" (i.e., all
of the data objects are contained within the root group). The file
contains:

+-----------------------+---------+------------------------------------+-----------------+
| Name                  | Type    | Contents                           | Attributess     |
+=======================+=========+====================================+=================+
| config                | group   | ShakeMap configuration             | None            |
+-----------------------+---------+------------------------------------+-----------------+
| info.json [1]_        | dataset | ShakeMap metadata in JSON format   | None            |
+-----------------------+---------+------------------------------------+-----------------+
| rupture.json [1]_     | dataset | GeoJSON representation of the      | None            |
|                       |         | finite fault                       |                 |
+-----------------------+---------+------------------------------------+-----------------+
| stationlist.json [1]_ | dataset | GeoJSON stationlist                | None            |
+-----------------------+---------+------------------------------------+-----------------+
| vs30                  | dataset | Vs30 values at output grid/points  | grid attrs [2]_ |
+-----------------------+---------+------------------------------------+-----------------+
| URATPGA [3]_          | dataset | Ratio of PGA stddev to GMPE stddev | grid attrs [2]_ |
+-----------------------+---------+------------------------------------+-----------------+
| IMT (multiple)        | dataset | Interpolated data for IMT          | grid attrs [2]_ |
+-----------------------+---------+------------------------------------+-----------------+
| IMT_sd (multiple)     | dataset | Stddev of interpolated IMT         | grid attrs [2]_ |
+-----------------------+---------+------------------------------------+-----------------+

.. [1] JSON datasets are stored as continuous strings.


.. [2] The atributes for grids are:

   +-----------+------------------------------------------------------+
   | Attr name | Contents                                             |
   +===========+======================================================+
   | type      | "grid"                                               |
   +-----------+------------------------------------------------------+
   | W         | The western boundary of the grid (degrees longitude) |
   +-----------+------------------------------------------------------+
   | E         | The eastern boundary of the grid (degrees longitude) |
   +-----------+------------------------------------------------------+
   | S         | The southern boundary of the grid (degrees latitude) |
   +-----------+------------------------------------------------------+
   | N         | The northern boundary of the grid (degrees latitude) |
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

.. [3] URATPGA is only provided for grids.


--------------------------------------

There will typically be multiple *IMT* (Intensity Measure Type) and 
*IMT_sd* (standard deviation of the corresponding IMT) datasets. For instance
*PGA*, *PGA_sd*, *PGV*, *PGV_sd*, *MMI*, *MMI_sd*, and various *SA(#num)* and
*SA(#num)_sd* [where #num is the period as a floating point number; e.g., 
*SA(1.0)*]. 

Developers will likely want to access **shake_result.hdf** through
the `shakelib OutputContainer class 
<https://usgs.github.io/shakelib/shakelib.utils.container.html>`_.
Also see, for example, *sm_contour* for some basic access patterns.
