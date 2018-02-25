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
the HDF5 data format can be found
`here <https://support.hdfgroup.org/HDF5/Tutor/HDF5Intro.pdf>`_.

*shake_result.hdf* consists of a number of groups and datasets. Our
implementation of HDF5 uses groups to contain Python dictionaries,
strings, and numpy arrays.  Dictionaries are stored as recursive groups.
For example, a Python *config* dictionary consisting of the following information:

.. code-block:: python

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

Each IMT dataset (MMI, PGA, etc.) is stored as a group containing two datasets, the mean values
for each cell and the standard deviations.  MMI data for the component 'Larger' will be stored
under a group called ``__imt_MMI_Larger__``. The mean array will be stored as
``mean``, and the standard deviation array will be stored as
``std``.  All IMT grid datasets will be accompanied by a dictionary of
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

All *shake_result.hdf* files will have a group ``__file_data_type__`` 
which will have a single attribute ``data_type`` that will be one of
'points' or 'grid'. This way the user can distinguish between the two
types of storage.

Regardless of whether the file stores grids or arrays of points, it will
also contain datasets of various distance parameters. These will be 
named ``distance_*`` where the wildcard will be replaced with  one of 
the typical 
source distance metrics (e.g., 'rrup' for rupture distance, 'rjb' for
Joyner-Boore distance, 'rhypo' for hypocentral distance, 'rx', 'ry0',
etc.) 
Similarly, the Vs30 data are found in a dataset named ``vs30``.
The metadata for distances and Vs30 consists of 'units' and 'digits'.

JSON datasets are stored as continuous strings.

There will typically be multiple *IMT* (Intensity Measure Type) datasets
(each containing the mean and standard deviation of the IMT). For instance
'PGA', 'PGV', 'MMI', and various 'SA(#num)' 
[where #num is the period as a floating point number; e.g., 
'SA(1.0)']. 

Python developers will likely want to access *shake_result.hdf* through
the `shakelib OutputContainer class 
<https://usgs.github.io/shakelib/shakelib.utils.container.html>`_.
Also see, for example, the *contour* module [:meth:`shakemap.coremods.contour`]
for some basic access patterns.

Matlab developers can use the function *read_shake_data.m*, which is included in
the repository for ShakeMap
`here <https://github.com/usgs/shakemap/blob/master/read_shakemap_data.m>`_.


Generic Amplification Factors
=============================

The ShakeMap generic amplification factor facility supports the inclusion
of linear amplifications that are not otherwise supported (by, for example,
Vs30-based
site amplifications), such as basin or topographic amplifications. The
ShakeMap operator may provide one or more files that contain factors
that will be added to the (natural logarithm) of the results returned
by the GMPE or IPE (the results from the IPE are not logged, but the
amplification factors are still additive). Mapped areas that extend
beyond the boundaries of the amplification factor file are given an
amplification factor of zero. If more than one amplification file is
present in the *GenericAmpFactors* directory, then the system will apply
all such files (i.e., the amplification factors will be cumulative to
the extent that the grids overlap).

The amplification factor file is a MapIO GridHDFContainer containing one
or more Grid2D objects corresponding to the IMTs to which they apply. For
instance, the following program creates a file **Test.hdf** which contains
grids for PGA, SA(0.3), SA(1.0), and SA(3.0). The grids are derived from 
GMT **.grd** files residing in the local directory::

    #! /usr/bin/env python

    from mapio.gmt import GMTGrid
    from mapio.gridcontainer import GridHDFContainer

    from shakelib.utils.imt_string import file_to_oq


    gc = GridHDFContainer.create('Test.hdf')

    files = ['PGA.grd', 'PSA0p3.grd', 'PSA1p0.grd', 'PSA3p0.grd']

    for myfile in files:
        g2d = GMTGrid.load(myfile)

        fbase, ext = myfile.split('.')
        name = file_to_oq(fbase)

        gc.setGrid(name, g2d)

    gc.close()

All of the grids in a given GridHDFContainer file must have exactly the same
boundaries and resolutions. The resulting HDF file should be placed in
*<install_dir>/data/GenericAmpFactors* where *<install_dir>* is the current
profile's install directory (as set/reported by **sm_profille**).

The rules for extracting and applying the amplification grids are as follows:

    - If an exact match to the output IMT is found, then that grid is used.
    - If the output IMT is 'SA(X)', where the period 'X' is between two of
      the SA periods in the amplifaction file, the grid that is applied 
      will the the 
      weighted average of the grids of the periods bracketing 'X'. The
      weighting will be the (normalized) log difference in the periods.
      I.e., if the bracketing periods are 'W' and 'Y", then the weight
      applied to the grid corresponding to period W ('gW') will be
      *wW = (log(Y) - log(X)) / (log(Y) - log(W))* and the weight for the
      grid corresponding to period Y ('gY') will be *wY = 1 - wW*, so
      the amplification factors used will be *wW * gW + wY * gY*.
    - If the period of the output IMT is less than the shortest period in
      the file, the grid corresponding to the shortest period will be used.
    - If the period of the output IMT is greater than the longest period
      in the file, the grid corresponding to the longest period will be used.
    - If the output IMT is PGA and PGA is not found in the file, it will be
      treated as SA(0.01) and the above rules will be applied.
    - If the output IMT is PGV and PGV is not found in the file, it will be
      treated as SA(1.0) and the above rules will be applied.
    - After the application of the above rules, if and IMT is not found, it will 
      be given amplification factors of zero.

Thus, if the output IMT is PGV, and PGV is not in the file, ShakeMap will
search for SA(1.0) using the rules above. If no SA grids are provided, the
resulting amplification grid will be all zeros.

If the operator wishes to alter these behaviors, then additional grids should
be included in the HDF file. For instance, if the extrapolation of the grids
for the longest and shortest periods to longer and shorter periods is 
undesirable, the operator should include grids (e.g., of zeros) just below 
and above the shortest and longest periods, respectively. If the interpolation
between periods is undesirable, then grids matching the output IMTs should be 
provided. Etc.


Stationlist GeoJSON
=============================

The *stationlist.json* file is a GeoJSON file describing the seismic station and
macroseismic data that comprised the input to the ShakeMap. In addition, the file
will contain predicted values and uncertainties for the station location from the 
selected GMPE, as well as the computed bias. The file also contains distance metrics,
and amplitudes converted from PGM to MMI or from MMI to PGM. 

To distinguish between seismic and macroseismic "stations", each station feature has,
within its **properties** section, and attribute **station_type**. The possible values are
**seismic** (for seismic instruments) and **macroseismic** (for "Did You Feel It?" or other
macroseismic observations.

The file consists of a list of "features," each representing one seismic station or
macroseismic observation. A typical seismic station with have a structure like this::

    {
      "id": "CI.CRG",
      "type": "Feature",
      "properties": {
        "intensity_flag": "",
        "network": "CI",
        "intensity_stddev": 0.63,
        "name": "Crocker Grade",
        "predictions": {
          "pgv": {
            "ln_phi": 0.603,
            "value": 0.1083,
            "ln_sigma": 0.7293,
            "ln_tau": 0.4101,
            "units": "cm/s",
            "ln_bias": -0.1386
          },
          "mmi": {
            "value": 2.181,
            "tau": 0.2618,
            "units": "intensity",
            "sigma": 0.7834,
            "phi": 0.7383,
            "bias": -0.1363
          },
          "pga": {
            "ln_phi": 0.6032,
            "value": 0.0363,
            "ln_sigma": 0.7401,
            "ln_tau": 0.4288,
            "units": "%g",
            "ln_bias": -0.7052
          },
          "sa(1.0)": {
            "ln_phi": 0.7328,
            "value": 0.1154,
            "ln_sigma": 0.8761,
            "ln_tau": 0.48,
            "units": "%g",
            "ln_bias": -0.0222
          },
          "sa(0.3)": {
            "ln_phi": 0.8036,
            "value": 0.0953,
            "ln_sigma": 1.0033,
            "ln_tau": 0.6007,
            "units": "%g",
            "ln_bias": -0.6593
          },
          "sa(3.0)": {
            "ln_phi": 0.6825,
            "value": 0.0319,
            "ln_sigma": 0.8201,
            "ln_tau": 0.4547,
            "units": "%g",
            "ln_bias": 0.4148
          }
        },
        "pga": 0.0106,
        "commType": "UNK",
        "mmi_from_pgm": {
          "sa(1.0)": {
            "value": 1,
            "sigma": 0.75
          },
          "pga": {
            "value": 1,
            "sigma": 0.66
          },
          "pgv": {
            "value": 1,
            "sigma": 0.63
          },
          "sa(3.0)": {
            "value": 1.58,
            "sigma": 0.89
          },
          "sa(0.3)": {
            "value": 1,
            "sigma": 0.82
          }
        },
        "pgv": 0.0174,
        "intensity": 1,
        "channels": [
          {
            "amplitudes": [
              {
                "value": 0.0063,
                "units": "%g",
                "name": "sa(3.0)",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0094,
                "units": "cm/s",
                "name": "pgv",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0056,
                "units": "%g",
                "name": "pga",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0108,
                "units": "%g",
                "name": "sa(0.3)",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0097,
                "units": "%g",
                "name": "sa(1.0)",
                "ln_sigma": 0,
                "flag": "0"
              }
            ],
            "name": "--.HHZ"
          },
          {
            "amplitudes": [
              {
                "value": 0.0057,
                "units": "%g",
                "name": "sa(3.0)",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0174,
                "units": "cm/s",
                "name": "pgv",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.01,
                "units": "%g",
                "name": "pga",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0244,
                "units": "%g",
                "name": "sa(0.3)",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0198,
                "units": "%g",
                "name": "sa(1.0)",
                "ln_sigma": 0,
                "flag": "0"
              }
            ],
            "name": "--.HHE"
          },
          {
            "amplitudes": [
              {
                "value": 0.0045,
                "units": "%g",
                "name": "sa(3.0)",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0108,
                "units": "cm/s",
                "name": "pgv",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0106,
                "units": "%g",
                "name": "pga",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0237,
                "units": "%g",
                "name": "sa(0.3)",
                "ln_sigma": 0,
                "flag": "0"
              },
              {
                "value": 0.0142,
                "units": "%g",
                "name": "sa(1.0)",
                "ln_sigma": 0,
                "flag": "0"
              }
            ],
            "name": "--.HHN"
          }
        ],
        "location": "",
        "distances": {
          "rx": 175.554,
          "rrup": 402.984,
          "rhypo": 403.17,
          "ry0": 363.216,
          "rjb": 403.183
        },
        "distance": 402.984,
        "source": "CI",
        "instrumentType": "OBSERVED",
        "station_type": "seismic",
        "code": "CI.CRG"
      },
      "geometry": {
        "coordinates": [
          -119.72486,
          35.2422
        ],
        "type": "Point"
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

A typical macroseismic "station" will have the following structure::

    {
      "properties": {
        "location": "",
        "distances": {
          "rrup": 24.255,
          "ry0": 22.916,
          "rx": -7.813,
          "rjb": 24.176,
          "rhypo": 36.231
        },
        "source": "DYFI",
        "channels": [
          {
            "amplitudes": [
              {
                "units": "intensity",
                "sigma": 0,
                "flag": "0",
                "name": "mmi",
                "value": 5.7
              }
            ],
            "name": "mmi"
          }
        ],
        "code": "294",
        "name": "UTM:(10S 0546 4261 1000)",
        "commType": "UNK",
        "instrumentType": "OBSERVED",
        "station_type": "macroseismic",
        "intensity_flag": "",
        "predictions": {
          "pgv": {
            "units": "cm/s",
            "ln_bias": -0.0068,
            "value": 9.3834,
            "ln_sigma": 0.6389,
            "ln_tau": 0.3458,
            "ln_phi": 0.5373
          },
          "sa(3.0)": {
            "units": "%g",
            "ln_bias": -0.0083,
            "value": 1.8265,
            "ln_sigma": 0.7397,
            "ln_tau": 0.3984,
            "ln_phi": 0.6233
          },
          "pga": {
            "units": "%g",
            "ln_bias": 0.0883,
            "value": 11.8117,
            "ln_sigma": 0.6354,
            "ln_tau": 0.3566,
            "ln_phi": 0.5259
          },
          "sa(0.3)": {
            "units": "%g",
            "ln_bias": -0.1608,
            "value": 29.1747,
            "ln_sigma": 0.6677,
            "ln_tau": 0.3267,
            "ln_phi": 0.5823
          },
          "mmi": {
            "units": "intensity",
            "sigma": 1.0797,
            "bias": -0.4454,
            "value": 5.7827,
            "tau": 0.4745,
            "phi": 0.9698
          },
          "sa(1.0)": {
            "units": "%g",
            "ln_bias": -0.0509,
            "value": 10.0253,
            "ln_sigma": 0.7533,
            "ln_tau": 0.3848,
            "ln_phi": 0.6476
          }
        },
        "intensity": 5.7,
        "distance": 24.255,
        "intensity_stddev": 0.3,
        "pga": 13.3845,
        "pgv": 8.9653,
        "pgm_from_mmi": {
          "pga": {
            "units": "%g",
            "value": 13.3845,
            "ln_sigma": 0.8059
          },
          "sa(0.3)": {
            "units": "%g",
            "value": 27.8831,
            "ln_sigma": 1.0131
          },
          "pgv": {
            "units": "cm/s",
            "value": 8.9653,
            "ln_sigma": 0.875
          },
          "sa(1.0)": {
            "units": "%g",
            "value": 9.1841,
            "ln_sigma": 1.0822
          },
          "sa(3.0)": {
            "units": "%g",
            "value": 1.9672,
            "ln_sigma": 1.4737
          }
        },
        "network": "DYFI"
      },
      "geometry": {
        "coordinates": [
          -122.4667,
          38.5007
        ],
        "type": "Point"
      },
      "type": "Feature",
      "id": "DYFI.294"
    }


The attributes of the macroseismic station are similar to those of the
seismic station (above), except:

- There will typically be only a single **channel** with a single **amplitude**
  element.
- The **pgm_from_mmi** section contains the output IMTs derived from MMI (to 
  the extent that the GMICE will make those conversions).
- Small intensity values (i.e., those less than 4.0) are not converted to
  PGM.

