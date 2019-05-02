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

+---------------+----------------------------+-------------+----------------------------------------------+
| Name          | Location                   | Python Type | Contents                                     |
+===============+============================+=============+==============================================+
| config        | /dictionaries/config       | JSON string | ShakeMap configuration                       |
+---------------+----------------------------+-------------+----------------------------------------------+
| info.json     | /dictionaries/info.json    | JSON string | ShakeMap metadata                            |
+---------------+----------------------------+-------------+----------------------------------------------+
| stations_dict | /dictionaries/stations_dict| JSON string | Dictionary representation of observed data   |
+---------------+----------------------------+-------------+----------------------------------------------+
| rupture       | /dictionaries/rupture      | JSON string | Dictionary representation of fault rupture   |
+---------------+----------------------------+-------------+----------------------------------------------+

It also will contain a number of arrays, which, when read with the HDFContainer getGrid() method, return
a Grid2D object which is a Python representation of a North-up 2D array of data, whose upper-left corner
coordinate and cell dimensions are known.  The definition of this object can be found here:
https://github.com/usgs/MapIO/blob/master/mapio/grid2d.py


Sampling of grids contained in the HDF:

+-------------------+----------------------------+-------------+----------------------------------------------+
| Name              | Location                   | Python Type | Contents                                     |
+===================+============================+=============+==============================================+
| vs30              | /arrays/vs30               | Grid2D      | Vs30 values at output grid/points            | 
+-------------------+----------------------------+-------------+----------------------------------------------+
| distance_rhypo    | /arrays/distances/rhypo    | Grid2D      | Hypocentral distance                         |
+-------------------+----------------------------+-------------+----------------------------------------------+
| distance_repi     | /arrays/distances/repi     | Grid2D      | Epicentral distance                          |
+-------------------+----------------------------+-------------+----------------------------------------------+
| distance_rjb      | /arrays/distances/rjb      | Grid2D      | RJB distance                                 |
+-------------------+----------------------------+-------------+----------------------------------------------+
| distance_rjb_std  | /arrays/distances/rjb_std  | Grid2D      | Standard deviation of the RJB distance if    |
|                   |                            |             | the point-source approximations are used,    |
|                   |                            |             | or will be zero if a finite fault is used.      |
+-------------------+----------------------------+-------------+----------------------------------------------+
| distance_rrup     | /arrays/distances/rrup     | Grid2D      | Rrup distance                                |
+-------------------+----------------------------+-------------+----------------------------------------------+
| distance_rrup_std | /arrays/distances/rrup_std | Grid2D      | Standard deviation of the Rrup distance if   |
|                   |                            |             | the point-source approximations are used,    |
|                   |                            |             | or will be zero if a finite fault is used.      |
+-------------------+----------------------------+-------------+----------------------------------------------+
| distance_rx       | /arrays/distances/rx       | Grid2D      | Rx distance (generalized coordinate used by  |
|                   |                            |             | some GMPEs)                                  |
+-------------------+----------------------------+-------------+----------------------------------------------+
| distance_ry0      | /arrays/distances/ry0      | Grid2D      | Ry0 distance (generalized coordinate used by |
|                   |                            |             | some GMPEs)                                  |
+-------------------+----------------------------+-------------+----------------------------------------------+


Each IMT dataset (MMI, PGA, etc.) is stored as a group containing two 
datasets: the mean values for each cell and the standard deviations.  
MMI data for the component 'Larger' will be stored under a group called 
``imts/MMI_GREATER_OF_TWO_HORIZONTAL``. The mean array will be stored as
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

+--------------+-----------------------------------------------------+-------------+---------------------+
| Name         | Location                                            | Python Type | Contents            |
+==============+=====================================================+=============+=============+=======+
| MMI Mean     | /arrays/imts/GREATER_OF_TWO_HORIZONTAL/MMI/mean     | Grid2D      | MMI Mean Values     | 
+--------------+-----------------------------------------------------+-------------+---------------------+
| MMI Std      | /arrays/imts/GREATER_OF_TWO_HORIZONTAL/MMI/std      | Grid2D      | MMI Std             | 
+--------------+-----------------------------------------------------+-------------+---------------------+
| Sa(0.3) Mean | /arrays/imts/GREATER_OF_TWO_HORIZONTAL/SA(0.3)/mean | Grid2D      | SA(0.3) Mean Values | 
+--------------+-----------------------------------------------------+-------------+---------------------+
| Sa(0.3) Std  | /arrays/imts/GREATER_OF_TWO_HORIZONTAL/SA(0.3)/std  | Grid2D      | SA(0.3) Std         | 
+--------------+-----------------------------------------------------+-------------+---------------------+

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

All *shake_result.hdf* files will have a group in the ``dictionaries``
group called ``file_data_type`` 
which will have a single attribute ``data_type`` that will be one of
'points' or 'grid'. This way the user can distinguish between the two
types of storage.

For grid files, there will also be sets of approximated attenuation
curves stored as one-dimensional arrays. These
will be labeled like ``/arrays/attenuation/<site>/<IMT>/<type>``
Where ``<IMT>`` will be one of the output IMTs (e.g., ``SA(3.0)``), 
``<site>`` will be one of ``rock`` or ``soil`` (for which ``rock``
means a Vs30 of 760 m/s, and ``soil`` means a Vs30 of 180 m/s), and
``<type>`` is one of ``mean`` (for the mean values) or ``std`` (for
the standard deviations). All units are in natural log space (except
for MMI). There will also be a set of datasets named like 
``/arrays/attenuation/distances/<type>`` which will contain the distances
(in km) corresponding to the points in the data arrays. The ``type`` will
be ``repi``, ``rhypo``, ``rjb``, ``rrup`` (for epicentral, hypocentral,
Joyner-Boore, and rupture distance, respectively). As with the other
distance arrays, ``rjb`` and ``rrup`` will be approximated if a finite
rupture model is not supplied in the input.

Dictionary datasets are stored as JSON strings.

There will typically be multiple *IMT* (Intensity Measure Type) datasets
(each containing the mean and standard deviation of the IMT). For instance
'PGA', 'PGV', 'MMI', and various 'SA(#num)' 
[where #num is the period as a floating point number; e.g., 
'SA(1.0)']. 

Python developers will likely want to access *shake_result.hdf* through
the ShakeMapOutputContainer class which may be found in the repository:
https://github.com/usgs/earthquake-impact-utils and imported from
``impactutils.utils.io.smcontainers``.
Also see, for example, the *contour* module [:meth:`shakemap.coremods.contour`]
for some basic access patterns.


.. _subsec-stationlist-geojson:

Stationlist GeoJSON
=============================

The *stationlist.json* file is a GeoJSON file describing the seismic station and
macroseismic data that comprised the input to the ShakeMap. It is contained within
*shake_result.hdf* and may be extracted into the *products* directory with the
**stations** module. In addition to the input station data, the file
will contain predicted values and uncertainties for the station location from the 
selected GMPE, as well as the computed bias. The file also contains distance metrics,
and amplitudes converted from PGM to MMI or from MMI to PGM. 

To distinguish between seismic and macroseismic "stations", each station feature has,
within its **properties** section, and attribute **station_type**. The possible values are
**seismic** (for seismic instruments) and **macroseismic** (for "Did You Feel It?" or other
macroseismic observations.

The file consists of a list of "features," each representing one seismic station or
macroseismic observation. A typical seismic station feature will have a structure 
like this::

    {
      "type": "Feature",
      "id": "NC.J051",
      "geometry": {
        "coordinates": [
          -122.007835,
          37.312901
        ],
        "type": "Point"
      },
      "properties": {
        "network": "NC",
        "intensity_flag": "",
        "mmi_from_pgm": [
          {
            "name": "sa(3.0)",
            "sigma": 0.89,
            "value": 3.75
          },
          {
            "name": "sa(1.0)",
            "sigma": 0.75,
            "value": 3.62
          },
          {
            "name": "sa(0.3)",
            "sigma": 0.82,
            "value": 3.19
          },
          {
            "name": "pgv",
            "sigma": 0.63,
            "value": 3.43
          },
          {
            "name": "pga",
            "sigma": 0.66,
            "value": 2.95
          }
        ],
        "distance": 104.211,
        "commType": "UNK",
        "intensity": 3.4,
        "pgv": 0.7679,
        "source": "NC",
        "instrumentType": "OBSERVED",
        "station_type": "seismic",
        "code": "NC.J051",
        "name": "So Tantau Av Cupertino",
        "pga": 0.4807,
        "intensity_stddev": 0.63,
        "distances": {
          "ry0": 103.951,
          "rrup": 104.211,
          "rjb": 104.208,
          "rx": 9.298,
          "rhypo": 104.433
        },
        "location": "",
        "channels": [
          {
            "amplitudes": [
              {
                "flag": "0",
                "units": "cm/s",
                "ln_sigma": 0,
                "name": "pgv",
                "value": 0.7679
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "sa(3.0)",
                "value": 0.2444
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "sa(1.0)",
                "value": 1.1346
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "pga",
                "value": 0.4807
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "sa(0.3)",
                "value": 1.1309
              }
            ],
            "name": "01.HNE"
          },
          {
            "amplitudes": [
              {
                "flag": "0",
                "units": "cm/s",
                "ln_sigma": 0,
                "name": "pgv",
                "value": 0.329
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "sa(3.0)",
                "value": 0.2168
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "sa(1.0)",
                "value": 0.5174
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "pga",
                "value": 0.2743
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "sa(0.3)",
                "value": 0.8392
              }
            ],
            "name": "01.HNZ"
          },
          {
            "amplitudes": [
              {
                "flag": "0",
                "units": "cm/s",
                "ln_sigma": 0,
                "name": "pgv",
                "value": 0.5312
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "sa(3.0)",
                "value": 0.2124
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "sa(1.0)",
                "value": 0.7154
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "pga",
                "value": 0.4429
              },
              {
                "flag": "0",
                "units": "%g",
                "ln_sigma": 0,
                "name": "sa(0.3)",
                "value": 1.1233
              }
            ],
            "name": "01.HNN"
          }
        ],
        "predictions": [
          {
            "units": "cm/s",
            "ln_sigma": 0.6356,
            "name": "pgv",
            "ln_phi": 0.5363,
            "value": 0.8747,
            "ln_bias": -0.1347,
            "ln_tau": 0.3412
          },
          {
            "units": "%g",
            "ln_sigma": 0.7032,
            "name": "pga",
            "ln_phi": 0.5689,
            "value": 1.186,
            "ln_bias": -0.7021,
            "ln_tau": 0.4134
          },
          {
            "units": "%g",
            "ln_sigma": 0.7337,
            "name": "sa(3.0)",
            "ln_phi": 0.6198,
            "value": 0.1489,
            "ln_bias": 0.4019,
            "ln_tau": 0.3927
          },
          {
            "units": "%g",
            "ln_sigma": 0.786,
            "name": "sa(0.3)",
            "ln_phi": 0.6556,
            "value": 2.3163,
            "ln_bias": -0.6296,
            "ln_tau": 0.4335
          },
          {
            "units": "%g",
            "ln_sigma": 0.7627,
            "name": "sa(1.0)",
            "ln_phi": 0.6539,
            "value": 0.7873,
            "ln_bias": -0.0214,
            "ln_tau": 0.3925
          },
          {
            "tau": 0.2178,
            "phi": 0.717,
            "units": "intensity",
            "bias": -0.1209,
            "name": "mmi",
            "value": 3.5145,
            "sigma": 0.7494
          }
        ]
      }
    }


The following features should be noted:

- The **coordinates** are given in longitude, latitude order.
- The units of the observed and predicted IMTs are provided; typically
  percent-g for accelerations and cm/s for velocity. The units of standard
  deviation and bias are in natural log units.
- **ln_tau** is the lagarithm of the between-even standard deviarion, **ln_phi**
  is the logarithm of the within-even standard deviation, and **ln_sigma**
  is the logarithm of the total standard deviation.
- Standard deviations for MMI are linear and omit the 'ln_' prefix.
- If the **flag** attribute is "0" or the empty string, the amplitude is 
  considered unflagged; any other value means the amplitude is flagged and
  therefore not included in the processing.
- The generic **distance** property is the same as **rrup** the rupture distance.
- The generic **intensity** property is the macroseismic intensity from the best
  available IMT.
- The **mmi_from_pgm** section contains the macroseismic intensity computed from
  the available IMTs (to the extent that the chosen GMICE is able to convert
  them).
- Floating point or integer values that cannot or were not determined will
  have the string value 'null'.

A typical macroseismic "station" feature will have the following structure::

    {
      "id": "DYFI.87",
      "type": "Feature",
      "geometry": {
        "type": "Point",
        "coordinates": [
          -122.6963,
          38.4474
        ]
      },
      "properties": {
        "intensity": 4.8,
        "predictions": [
          {
            "units": "intensity",
            "name": "mmi",
            "sigma": 1.0851,
            "value": 5.1036,
            "phi": 0.9733,
            "tau": 0.4796,
            "bias": -0.4463
          },
          {
            "name": "sa(0.3)",
            "ln_bias": -0.1675,
            "value": 18.2415,
            "ln_sigma": 0.7003,
            "ln_tau": 0.3563,
            "ln_phi": 0.6029,
            "units": "%g"
          },
          {
            "name": "sa(1.0)",
            "ln_bias": -0.0512,
            "value": 6.0597,
            "ln_sigma": 0.7585,
            "ln_tau": 0.389,
            "ln_phi": 0.6511,
            "units": "%g"
          },
          {
            "name": "sa(3.0)",
            "ln_bias": -0.0083,
            "value": 1.0917,
            "ln_sigma": 0.7376,
            "ln_tau": 0.3964,
            "ln_phi": 0.622,
            "units": "%g"
          },
          {
            "name": "pgv",
            "ln_bias": -0.0068,
            "value": 5.721,
            "ln_sigma": 0.6437,
            "ln_tau": 0.3495,
            "ln_phi": 0.5406,
            "units": "cm/s"
          },
          {
            "name": "pga",
            "ln_bias": 0.0897,
            "value": 7.5028,
            "ln_sigma": 0.6602,
            "ln_tau": 0.3775,
            "ln_phi": 0.5416,
            "units": "%g"
          }
        ],
        "distance": 35.27,
        "pgv": 4.5832,
        "pga": 6.8063,
        "pgm_from_mmi": [
          {
            "value": 1.0441,
            "ln_sigma": 1.4737,
            "name": "sa(3.0)",
            "units": "%g"
          },
          {
            "value": 4.7097,
            "ln_sigma": 1.0822,
            "name": "sa(1.0)",
            "units": "%g"
          },
          {
            "value": 4.5832,
            "ln_sigma": 0.875,
            "name": "pgv",
            "units": "cm/s"
          },
          {
            "value": 6.8063,
            "ln_sigma": 0.8059,
            "name": "pga",
            "units": "%g"
          },
          {
            "value": 14.9458,
            "ln_sigma": 1.0131,
            "name": "sa(0.3)",
            "units": "%g"
          }
        ],
        "channels": [
          {
            "amplitudes": [
              {
                "value": 4.8,
                "name": "mmi",
                "flag": "0",
                "sigma": 0,
                "units": "intensity"
              }
            ],
            "name": "mmi"
          }
        ],
        "intensity_stddev": 0.3,
        "name": "UTM:(10S 0526 4255 1000)",
        "instrumentType": "OBSERVED",
        "commType": "UNK",
        "location": "",
        "distances": {
          "rrup": 35.27,
          "ry0": 20.571,
          "rjb": 35.219,
          "rx": -28.528,
          "rhypo": 43.728
        },
        "network": "DYFI",
        "intensity_flag": "",
        "station_type": "macroseismic",
        "code": "87",
        "source": "DYFI"
      }
    }

The attributes of the macroseismic station are similar to those of the
seismic station (above), except:

- There will typically be only a single **channel** with a single **amplitude**
  element.
- The **pgm_from_mmi** section contains the output IMTs derived from MMI (to 
  the extent that the GMICE will make those conversions).
- Small intensity values (i.e., those less than 4.0) are not converted to
  PGM (i.e., they will have the value 'null').

