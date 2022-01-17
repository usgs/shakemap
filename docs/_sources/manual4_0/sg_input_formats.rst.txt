.. _sec-input-formats-4:

****************************
Input Data Formats
****************************

ShakeMap has a number of input files, some of which may be provided in more 
than one format. See below for a discussion of some of them.

.. _subsec-event-xml-file:

Event XML File
=============================

This file is *required*.

The earthquake origin information is contained in a simple XML file
called *event.xml* that is contained in the event's *current*
directory. An example file::

  <earthquake id="us1000db5t" netid="us" lat="38.7161" lon="69.9779"
  depth="5.0" mag="5.7" locstring="11km SE of Roghun, Tajikistan"
  network="USGS National Earthquake Information Center, PDE"
  time="2018-03-29T22:54:12Z" mech="SS" reference="Smith, et. al. 2018"/>

The attributes of the earthquake element are:

+-----------------------+-------------------------------------------------------+
| Field                 | Description                                           |
+=======================+=======================================================+
| id                    | A string, consisting of the event ID. This should be  |
|                       | be unique for the network code **netid** (see below). |
+-----------------------+-------------------------------------------------------+
| netid                 | Usually a two character network code, but can be any  |
|                       | string.                                               |
+-----------------------+-------------------------------------------------------+
| network               | Long description of the organization behind netid.    |
|                       | This field is required, but may be the empty string   |
+-----------------------+-------------------------------------------------------+
| lat                   | Latitude at earthquake hypocenter.                    |
+-----------------------+-------------------------------------------------------+
| lon                   | Longitude at earthquake hypocenter.                   |
+-----------------------+-------------------------------------------------------+
| depth                 | Depth at earthquake hypocenter.                       |
+-----------------------+-------------------------------------------------------+
| mag                   | Magnitude of earthquake.                              |
+-----------------------+-------------------------------------------------------+
| time                  | Time of earthquake origin, always UTC in              |
|                       | YYYY-MM-DDTHH:MM:SS.fZ format.                        |
+-----------------------+-------------------------------------------------------+
| locstring             | String describing earthquake location.                |
+-----------------------+-------------------------------------------------------+
| mech                  | (optional) Focal mechanism, one of "RS", "SS",        |
|                       | "NM", or "ALL".                                       |
+-----------------------+-------------------------------------------------------+
| reference             | (optional) Source for hypocenter information          |
+-----------------------+-------------------------------------------------------+
| event_type            | (optional) Event type, one of 'ACTUAL' or 'SCENARIO'  |
+-----------------------+-------------------------------------------------------+
| productcode           | (optional) Used to distinguish between different      |
|                       | ShakeMaps created for the same event (e.g., one map   |
|                       | showing the extent of shaking, and another zoomed     |
|                       | into a city of interest). If the productcode is not   |
|                       | supplied, it will be taken from the name of the       |
|                       | data directory (which is also the "event_id" that is  |
|                       | supplied to the **shake** program).                   |
+-----------------------+-------------------------------------------------------+

Moment XML File
=============================

This file is *optional*.

Users may provide moment tensor information in a QuakeML file called
*moment.xml* that is in the event's *current* directory. The QuakeML
format is described exhaustively here: https://quake.ethz.ch/quakeml/

The QuakeML *must* contain a *focalMechanism* element with at least
complete *nodalPlanes* and *principalAxes* elements. See the QuakeML
documentation for details.

ShakeMap XML Data File
======================

These files are *optional*.

To inform your ShakeMap with real data, either macroseismic intensity
or instrumented accelerations and velocities, you may include any
number of XML files ending in "\*_dat.xml" as well as one called
"stationlist.xml" in the event's *current*
directory. These files consists of a number of *elements*, with
defined *attributes*.

The *earthquake* element (optional) is described above in the section
:ref:`subsec-event-xml-file`.

Following the earthquake element is the *stationlist* element, which
has a *created* attribute, which is the same as the one in the
earthquake element.

Example stationlist element::

   <stationlist created="1529614595">

The stationlist element contains zero to many *station* elements, with the following attributes:

+-------------+-----------------------------------------------------------------+
| Attribute   | Description                                                     |
+=============+=================================================================+
| code        | The station code.                                               |
+-------------+-----------------------------------------------------------------+
| name        | The station name and/or description.                            |
+-------------+-----------------------------------------------------------------+
| insttype    | Description of instrument type.                                 |
+-------------+-----------------------------------------------------------------+
| lat         | Station latitude (in decimal degrees).                          |
+-------------+-----------------------------------------------------------------+
| lon         | Station longitude (in decimal degrees).                         |
+-------------+-----------------------------------------------------------------+
| source      | Agency that maintains the station (i.e., SCSN, NSMP ,...).      |
+-------------+-----------------------------------------------------------------+
| netid       | The network ID string; for MMI observations this must be one    |
|             | of ‘MMI,’ ‘CIIM,’ ‘DYFI,’ or ‘INTENSITY.’.                      |
+-------------+-----------------------------------------------------------------+
| commtype    | Digital or analog communications (DIG or ANA).                  |
+-------------+-----------------------------------------------------------------+
| loc         | (optional) Free form text describing station location.          |
+-------------+-----------------------------------------------------------------+
| intensity   | The intensity value of the observation (decimal).               |
+-------------+-----------------------------------------------------------------+

Example station element::

   <station code="ADO" name="Adelanto Receiving Station"
   insttype="TriNet" lat="34.55046" lon="-117.43391" source="SCSN
   and TriNet" commtype="DIG" netid=”CI” loc="Adelanto, on Hwy 395">

The station element can contain zero to many *comp* elements. Each
comp element can contain one *acc* element and one *vel* element,
and may contain *psa03*, *psa10* and *psa30* elements (one of
each). These refer to peak acceleration, velocity, and 5%-damped
pseduo-spectral acceleration (at 0.3, 1.0, and 3.0 sec period) values
for the named channel at the named station. The acc, vel, psa03,
psa10, and psa30 elements are empty but have the following attributes:

+-------------+-----------------------------------------------------------------+
| Attribute   | Description                                                     |
+=============+=================================================================+
| value       | The amplitude value.                                            |
+-------------+-----------------------------------------------------------------+
| flag        | Flag indicating problematic data (optional).                    |
+-------------+-----------------------------------------------------------------+
| units       | (optional) See below.                                           |
+-------------+-----------------------------------------------------------------+

If the units are unspecifies, the value attributes are expected to have units of:

+-------------+-----------------------------------------------------------------+
| Attribute   | Units                                                           |
+=============+=================================================================+
| acc, pga    | %g (i.e., percent of the earth’s nominal gravitational          |
|             | acceleration).                                                  |
+-------------+-----------------------------------------------------------------+
| vel, pgv    | cm/s (centimeters per second).                                  |
+-------------+-----------------------------------------------------------------+
| psa         | %g.                                                             |
+-------------+-----------------------------------------------------------------+

The **units** attribute may also be specified if the amplitudes are in
logarithmic (natural) units:

+-------------+-----------------------------------------------------------------+
| Attribute   | Units Designator                                                |
+=============+=================================================================+
| acc, pga    | "ln(g)" (i.e., log of the earth’s nominal gravitational         |
|             | acceleration).                                                  |
+-------------+-----------------------------------------------------------------+
| vel, pgv    | "ln(cm/s)" (log of centimeters per second).                     |
+-------------+-----------------------------------------------------------------+
| psa         | "ln(g)".                                                        |
+-------------+-----------------------------------------------------------------+

The operator may also specify a standard deviations for "observations" that are
the mean of a distribution. This standard deviation is specified in natural
logarithmic units regardless of the units of the amplitudes themselves. The
standard deviations are specified with the **ln_sigma** attribute.

The flag attribute indicates problematic data. Any value other than
“0” (zero) or “” (i.e., an empty string) will cause ShakeMap to reject
the amplitude (and, in fact, all the amplitudes of that type for that
station). Though any non-zero flag will kill an amplitude, the
following flags are currently defined:

+-------------+-----------------------------------------------------------------+
| Flag        | Description                                                     |
+=============+=================================================================+
| T           | Automatically flagged by ShakeMap as an outlier.                |
+-------------+-----------------------------------------------------------------+
| M           | Manually flagged (in grind.conf) by the ShakeMap operator.      |
+-------------+-----------------------------------------------------------------+
| G           | Glitch. Amplitude clipped or below instrument noise threshold.  |
+-------------+-----------------------------------------------------------------+
| I           | Incomplete (a data gap existed in the time window used to       |
|             | calculate the amplitude).                                       |
+-------------+-----------------------------------------------------------------+

An abbreviated example of a complete station data file::

  <?xml version="1.0" encoding="UTF-8" standalone="yes"?>
  <!DOCTYPE stationlist [
  ... DTD description ...
  ]>
  <stationlist created="1070030689">
  <station code="ADO" name="Adelanto Receiving Station"
  insttype="TriNet" lat="34.55046" lon="-117.43391" source="SCSN
  and TriNet" commtype="DIG" netid=”CI” loc="Adelanto, on Hwy 395
  ">
  <comp name="HHE">
  <acc value="0.0083" flag="0" />
  <vel value="0.0030" flag="0" />
  <psa03 value="0.0146" flag="0" />
  <psa10 value="0.0049" flag="0" />
  <psa30 value="0.0003" flag="0" />
  </comp>
  <comp name="HHN">
  <acc value="0.0088" flag="0" />
  <vel value="0.0028" flag="0" />
  <psa03 value="0.0111" flag="0" />
  <psa10 value="0.0040" flag="0" />
  <psa30 value="0.0004" flag="0" />
  </comp>
  <comp name="HHZ">
  <acc value="0.0087" flag="0" />

  <vel value="0.0016" flag="0" />
  <psa03 value="0.0080" flag="0" />
  <psa10 value="0.0013" flag="0" />
  <psa30 value="0.0002" flag="0" />
  </comp>
  </station>
  ... additional station tags (omitted)...
  <station code="WSS" name="West Side Station" insttype="TriNet"
  lat="34.1717" lon="-118.64971" source="SCSN and TriNet"
  commtype="DIG" netid=”CI” loc="Hidden Hills, Valley Circle Dr.">
  <comp name="HHE">
  <acc value="0.0225" flag="0" />
  <vel value="0.0031" flag="0" />
  <psa03 value="0.0182" flag="0" />
  <psa10 value="0.0016" flag="0" />
  <psa30 value="0.0002" flag="0" />
  </comp>
  <comp name="HHN">
  <acc value="0.0209" flag="0" />
  <vel value="0.0029" flag="0" />
  <psa03 value="0.0234" flag="0" />
  <psa10 value="0.0019" flag="0" />
  <psa30 value="0.0001" flag="0" />
  </comp>
  <comp name="HHZ">
  <acc value="0.0187" flag="0" />
  <vel value="0.0020" flag="0" />
  <psa03 value="0.0073" flag="0" />
  <psa10 value="0.0005" flag="0" />
  <psa30 value="0.0000" flag="0" />
  </comp>
  </station>
  </stationlist>

Intensity data uses the same format of input XML as other ground
motion data, but uses three new attributes to the station tag: the
**intensity** attribute should be set to the decimal intensity for the
“station;” the **intensity_stddev** should specify the standard deviation
of the intensity observation; the **intensity_flag** should specify the 
flag (usually "0") of the observation (see flag table, above). Also the
netid attribute should be set to “MMI,” “CIIM,” “DYFI,”
or “INTENSITY” (all four are currently equivalent). If netid is set to
one of these values, any amplitude data (i.e., data enclosed in a comp
tag) will be ignored and *model* will use the configured GMICE to derive
the ground motions. Likewise, if netid is not one of these values, the
intensity attribute will be ignored and grind will compute intensity
using the GMICE.

Below is an example of a station tag that contains intensity information::

  <station code="91042" name="ZIP Code 91042 (Intensity VII, 38
  responses)" insttype="USGS (Did You Feel It?)" lat="34.282604"
  lon="-118.237943" source="USGS (Did You Feel It?)" netid="CIIM"
  commtype="USGS (Did You Feel It?)" intensity="7.4" intensity_stddev="0.3"
  intensity_flag="0">

The earthquake and stationlist XML files are combined in the GeoJSON
output file provided to the public. 

.. _subsec-json-input-stations-4:

ShakeMap JSON Data File
=======================

ShakeMap will also accept a ShakeMap-produced GeoJSON *stationlist.json*
file as input (see :ref:`subsec-stationlist-geojson`). Additional 
JSON files of the form *\*_dat.json* file may also be included in the input.

The information contained in the JSON input files is similar to that in
the XML input files (see above), but is structured differently::

    {
      "type": "FeatureCollection",
      "features": [
        {
          "geometry": {
            "type": "Point",
            "coordinates": [
              143.157196,
              42.014999
            ]
          },
          "type": "Feature",
          "id": "II.ERM",
          "properties": {
            "name": "Erimo, Hidaka, Hokkaido, Japan",
            "code": "II.ERM",
            "pgv": "null",
            "commType": "UNK",
            "vs30": 760,
            "intensity": "null",
            "network": "II",
            "distance": 462.284,
            "source": "II",
            "channels": [
              {
                "amplitudes": [
                  {
                    "name": "sa(3.0)",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0009,
                    "units": "%g"
                  },
                  {
                    "name": "pgv",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0056,
                    "units": "cm/s"
                  },
                  {
                    "name": "sa(1.0)",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0051,
                    "units": "%g"
                  },
                  {
                    "name": "pga",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0118,
                    "units": "%g"
                  },
                  {
                    "name": "sa(0.3)",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0201,
                    "units": "%g"
                  }
                ],
                "name": "BHZ"
              },
              {
                "amplitudes": [
                  {
                    "name": "sa(3.0)",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.001,
                    "units": "%g"
                  },
                  {
                    "name": "pgv",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0058,
                    "units": "cm/s"
                  },
                  {
                    "name": "sa(1.0)",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0069,
                    "units": "%g"
                  },
                  {
                    "name": "pga",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0146,
                    "units": "%g"
                  },
                  {
                    "name": "sa(0.3)",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.026,
                    "units": "%g"
                  }
                ],
                "name": "BH2"
              },
              {
                "amplitudes": [
                  {
                    "name": "sa(3.0)",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0012,
                    "units": "%g"
                  },
                  {
                    "name": "pgv",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0073,
                    "units": "cm/s"
                  },
                  {
                    "name": "sa(1.0)",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0046,
                    "units": "%g"
                  },
                  {
                    "name": "pga",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0182,
                    "units": "%g"
                  },
                  {
                    "name": "sa(0.3)",
                    "ln_sigma": 0,
                    "flag": "0",
                    "value": 0.0235,
                    "units": "%g"
                  }
                ],
                "name": "BH1"
              }
            ],
            "station_type": "seismic",
            "intensity_flag": "",
            "location": "",
            "intensity_stddev": "null",
            "instrumentType": "OBSERVED",
          }
        },
        <additional "features" (i.e., stations)>
      ]
    }

Note that the names of the intensity measure types are lower case,
and the spectral accelerations are of the form *sa(1.0)* where the
number in paraentheses is the period. Additional fields may be present
in the JSON file, but they will be ignored. Intensity observations
should have a **netid** as specified for the XML files (see above),
and should have a **channels** element that is an empty list
(i.e., "channels: []").

Source Text File
================

Because most ShakeMap installations automatically generate XML input
files and write them to the input directory, manual changes made by
the operator to the event.xml file will generally be overwritten by
the next automatic run. We therefore provide a mechanism by which the
operator may override or supplement any of the event-specific data in
event.xml. The operator may add an optional file to an event’s input
directory called *source.txt*. The structure of the file is one
parameter per line, in the form *parameter=value*. In particular, the
operator may specify the source mechanism with “mech” (this is the
equivalent of the “type” attribute in event.xml), which may be one of
“RS,” “SS,” “NM,” or “ALL” for reverse slip, strike slip, normal, and
unspecified mechanisms, respectively.  Any of the other source
parameters may also be set: eid, location, time, lat, lon, depth, mag,
etc.. Blank lines and lines beginning with ‘#’ (i.e., comments) are
ignored.


Rupture Specification
=====================
There are three classes of rupture objects:

- `PointRupture`
- `QuadRupture`
- `EdgeRupture`
  
A `PointRupture` is just a point representation of the earthquake and
is generated from the origin and so no additional specification is
required. In this case, distance calculations use approximate adjustments
to convert from epicentral distance to finite distances based on the
earthquake magnitude.

There are two extended-source rupture objects: a `QuadRupture` and an
`EdgeRupture`. There are some general rules for the specification of the
rupture vertices that apply to both of the extended-source rupture
objects:

- Vertices must start on the top edge of the rupture.
- The top and bottom edges must contain the same number of vertices.
- The first and last points must be identical to close the polygon, and this
  means that there must always be an odd number of vertices.
- The top edge of the rupture must always be above the bottom edge.


In cross section, a single-segment multiple-quadrilateral rupture might look
schematically like this::

       _.-P1-._
    P0'        'P2---P3
    |                  \
    P7---P6----P5-------P4

An `EdgeRupture` does not have any additional constraints beyond those already
described. This rupture would be initialized as an `EdgeRupture` because P1 causes
the top edges of two of the constituent quadrilaterals to not be horizontal, which
is a requirement for the `QuadRupture` class.

A `QuadRupture` consists of one or more quadrilaterals that can be grouped
into segments. The distance calculations are faster for a `QuadRupture` than
and `EdgeRupture`. The additional requirements for the vertices of a `QuadRupture`
are:

- The top and bottom edges of each quadrilateral are horizontal. In the example
  there are three quadrilateriasl: P0-P1-P6-P7, P1-P2-P5-P6, P2-P3-P4-P5. Of those,
  only the last one fulfills this criteria.
- The four points that define each quadrilateral must be approximately co-planar.

In ShakeMap version 3, ruptures were specified in a `*_fault.txt` file format. We
still support this format for backwards compatibility but we prefer to use the
GeoJSON format described below. Eventually we will stop support for the older
`*_fault.txt` file format.

Rupture GeoJson File
====================

This file is *optional*.

Rupture (also referred to as "finite fault") files are defined in ShakeMap 4
as GeoJSON files, a standard format for representing geospatial
data. This format is described in great detail here:
https://tools.ietf.org/html/rfc7946

The rupture format consists of a *FeatureCollection*, containing one
to many Features. The FeatureCollection should contain a dictionary
called *metadata*, which contains of the following fields:

+-----------------------+-------------------------------------------------------+
| Field                 | Description                                           |
+=======================+=======================================================+
| id                    | A unique string, consisting of netid plus             |
|                       | event ID.                                             |
+-----------------------+-------------------------------------------------------+
| netid                 | Usually a two character network code, but can be any  |
|                       | string.                                               |
+-----------------------+-------------------------------------------------------+
| network               | Long description of the organization behind netid.    |
+-----------------------+-------------------------------------------------------+
| lat                   | Latitude at earthquake hypocenter.                    |
+-----------------------+-------------------------------------------------------+
| lon                   | Longitude at earthquake hypocenter.                   |
+-----------------------+-------------------------------------------------------+
| depth                 | Depth at earthquake hypocenter.                       |
+-----------------------+-------------------------------------------------------+
| mag                   | Magnitude of earthquake.                              |
+-----------------------+-------------------------------------------------------+
| time                  | Time of earthquake origin, always UTC in              |
|                       | YYYY-MM-DDTHH:MM:SSZ format.                          |
+-----------------------+-------------------------------------------------------+
| locstring             | String describing earthquake location.                |
+-----------------------+-------------------------------------------------------+
| reference             | Source for rupture information.                       |
+-----------------------+-------------------------------------------------------+
| mech                  | Focal mechanism, one of "RS", "SS",                   |
|                       | "NM", or "ALL".                                       |
+-----------------------+-------------------------------------------------------+

Note that the only *required* field for specifying a rupture is *reference* and
that the other fields are merged with origin information and are included when
this file is output after running ShakeMap.

Each Feature must contain either a *Point* or *MultiPolygon*
geometry. Note that there is usually no reason to use the *Point* Feature type
when specifying a rupture, but the output rupture file is a *Point* type for
`PointRupture` object. 


The file should be named *rupture.json* and placed in the event's
*current* directory. Here is a single-segment single-quadrilaterial example::

  {
    "type": "FeatureCollection",
    "metadata": {
      "reference": "Wald, D. J., T. H. Heaton, and K. W. Hudnut (1996). The Slip History of the 1994 Northridge, California, Earthquake Determined from Strong-Motion, Teleseismic, GPS, and Leveling Data, Bull. Seism. Soc. Am. 86, S49-S70."
    },
    "features": [{
        "type": "Feature",
        "properties": {
          "rupture type": "rupture extent"
        },
        "geometry": {
          "type": "MultiPolygon",
          "coordinates": [
            [
              [
                [-118.421, 34.315, 5.0],
                [-118.587, 34.401, 5.0],
                [-118.693, 34.261, 20.427],
                [-118.527, 34.175, 20.427],
                [-118.421, 34.315, 5.0]
              ]
            ]
          ]
        }
    }]
  }

Here is a single-segment multi-quadrilaterial example::

    {
      "type": "FeatureCollection",
      "metadata": {
        "reference": "Konca, A. O, Hjorleifsdottir, V., Song, T. A., Avouac, J., Helmberger, D., Ji, C., Sieh, K., Briggs, R., and A. Meltzner. Rupture Kinematics of the 2005 Mw 8.6 Nias-Simeulue Earthquake from the Joint Inversion of Seismic and Geodetic Data (2007). BSSA Vol. 97, No. 1A, pp. S307-S322, January 2007, doi: 10.1785/0120050632."
      },
      "features": [{
          "type": "Feature",
          "properties": {
            "rupture type": "rupture extent"
          },
          "geometry": {
            "type": "MultiPolygon",
            "coordinates": [
              [
                [
                  [97.8322, -0.0442476, 10],
                  [96.5212, 0.897694, 10],
                  [95.9014, 2.57177, 10],
                  [96.9539, 3.18268, 40],
                  [98.1601, 1.89031, 40],
                  [98.7634, 1.09946, 40],
                  [97.8322, -0.0442476, 10]
		]
              ]
            ]
          }
      }]
    }

Here is a multi-segment example::

    {
      "type": "FeatureCollection",
      "metadata": {
        "reference": "Oglesby, D. D., D. S. Dreger, R. A. Harris, N. Ratchkovski, and R. Hansen (2004). Inverse kinematic and forward dynamic models of the 2002 Denali fault earthquake, Alaska, Bull. Seism. Soc. Am. 94, S214-S233."
      },
      "features": [{
          "type": "Feature",
          "properties": {
            "rupture type": "rupture extent"
          },
          "geometry": {
            "type": "MultiPolygon",
            "coordinates":[
              [
                [
                  [-147.807, 63.434, 0.0],
                  [-147.21, 63.472, 0.0],
                  [-147.267, 63.65, 22.294],
                  [-147.864, 63.613, 22.294],
                  [-147.807, 63.434, 0.0]
                ],
                [
                  [-146.951, 63.551, 0.0],
                  [-147.551, 63.518, 0.0],
                  [-147.551, 63.518, 30.0],
                  [-146.951, 63.551, 30.0],
                  [-146.951, 63.551, 0.0]
                ],
                [
                  [-145.968, 63.453, 0.0],
                  [-146.952, 63.547, 0.0],
                  [-146.952, 63.547, 30.0],
                  [-145.968, 63.453, 30.0],
                  [-145.968, 63.453, 0.0]
                ],
                [
                  [-143.586, 62.872, 0.0],
                  [-145.996, 63.427, 0.0],
                  [-145.996, 63.427, 30.0],
                  [-143.586, 62.872, 30.0],
                  [-143.586, 62.872, 0.0]
                ],
                [
                  [-142.5, 62.114, 0.0],
                  [-143.669, 62.831, 0.0],
                  [-143.669, 62.831, 30.0],
                  [-142.5, 62.114, 30.0],
                  [-142.5, 62.114, 0.0]
                ]
              ]
            ]
          }
      }]
    }

    
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
      the SA periods in the amplification file, the grid that is applied 
      will be the 
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

