.. _sec-input-formats-4:

****************************
Input Data Formats
****************************

Event XML File
=============================

This file is *required*.

The earthquake origin information is contained in a simple XML file
called *event.xml* that is contained in the event's *current*
directory. An example file::

  <earthquake id="us1000db5t" netid="us" network="USGS National Earthquake Information Center, PDE"
  lat="38.7161" lon="69.9779" depth="5.0" mag="5.7" time="2018-03-29T22:54:12Z"
  locstring="11km SE of Roghun, Tajikistan" mech="SS" reference="Smith, et. al. 2018"/>

The attributes of the earthquake element are:

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
| mech                  | (optional) Focal mechanism, one of "RS", "SS",        |
|                       | "NM", or "ALL".                                       |
+-----------------------+-------------------------------------------------------+
| reference             | (optional) Source for hypocenter information          |
+-----------------------+-------------------------------------------------------+
| productcode           | (optional) Used to distinguish between different      |
|                       | ShakeMaps created for the same event (i.e., one map   |
|                       | showing the extent of shaking, and another zoomed     |
|                       | into a city of interest.                              |
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
number of XML files ending in "*_dat.xml" in the event's *current*
directory. These files consists of a number of *elements*, with
defined *attributes*.

The *earthquake* element (optional) has the following attributes:

+-------------+-----------------------------------------------------------------+
| Attribute   | Description                                                     |
+=============+=================================================================+
| id          | The event id.                                                   |
+-------------+-----------------------------------------------------------------+
| created     | File creation time (Unix epoch -- seconds since Jan 1, 1970).   |
+-------------+-----------------------------------------------------------------+
| lat         | Latitude (in decimal degrees, negative in southern hemisphere). |
+-------------+-----------------------------------------------------------------+
| lon         | Longitude (in decimal degrees, negative in western hemisphere). |
+-------------+-----------------------------------------------------------------+
| depth       | Depth in km, positive down                                      |
+-------------+-----------------------------------------------------------------+
| locstring   | a free-form descriptive string of location relative to          |
|             | landmarks                                                       |
+-------------+-----------------------------------------------------------------+
| mag         | Earthquake magnitude.                                           |
+-------------+-----------------------------------------------------------------+
| type        | A string specifying the rupture type; the accepted types are    |
|             | RS, SS, NM, and ALL, for reverse slip, strike slip, normal,     |
|             | and unspecified ruptures, respectively.                         |
+-------------+-----------------------------------------------------------------+
| year        | Year, 4 digit format.                                           |
+-------------+-----------------------------------------------------------------+
| month       | Month, 1-12.                                                    |
+-------------+-----------------------------------------------------------------+
| day         | Day, 1-31.                                                      |
+-------------+-----------------------------------------------------------------+
| hour        | Hour, 0-23.                                                     |
+-------------+-----------------------------------------------------------------+
| minute      | Minute, 0-59.                                                   |
+-------------+-----------------------------------------------------------------+
| second      | Second, 0-59.                                                   |
+-------------+-----------------------------------------------------------------+
| timezone    | Timezone abbreviation (i.e., GMT, PST, PDT).                    |
+-------------+-----------------------------------------------------------------+

Example earthquake element::
  
   <earthquake id="14000376" lat="34.2722" lon="-118.7530"
   mag="3.6" year="2003" month="10" day="29" hour="23" minute="44"
   second="48" timezone="GMT" depth="13.81" locstring="2.6 mi W of
   Simi Valley, CA" created="1069292035" />

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
comp element must contain one *acc* element, and one *vel* element,
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

The value attributes are expected to have units of:

+-------------+-----------------------------------------------------------------+
| Attribute   | Units                                                           |
+=============+=================================================================+
| acc         | %g (i.e., percent of the earth’s nominal gravitational          |
|             | acceleration).                                                  |
+-------------+-----------------------------------------------------------------+
| vel         | cm/s (centimeters per second).                                  |
+-------------+-----------------------------------------------------------------+
| psa         | %g.                                                             |
+-------------+-----------------------------------------------------------------+

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
| G           | Amplitude clipped or below the instrument noise threshold.      |
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
motion data, but uses two new attributes to the station tag: the
intensity attribute should be set to the decimal intensity for the
“station;” the netid attribute should be set to “MMI,” “CIIM,” “DYFI,”
or “INTENSITY” (all four are currently equivalent). If netid is set to
one of these values, any amplitude data (i.e., data enclosed in a comp
tag) will be ignored and grind will use the mmi2pgm function to derive
the ground motions. Likewise, if netid is not one of these values, the
intensity attribute will be ignored and grind will compute intensity
using the pgm2mmi function.

Below is an example of a station tag that contains intensity information::

  <station code="91042" name="ZIP Code 91042 (Intensity VII, 38
  responses)" insttype="USGS (Did You Feel It?)" lat="34.282604"
  lon="-118.237943" source="USGS (Did You Feel It?)" netid="CIIM"
  commtype="USGS (Did You Feel It?)" intensity="7.4">

The earthquake and stationlist XML files are combined in the GeoJSON
output file provided to the public. 


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

Rupture GeoJson File
====================

This file is *optional*.

Rupture (also known as finite fault) files are defined in ShakeMap 4
as GeoJSON files, a standard format for representing geospatial
data. This format is described in great detail here: https://tools.ietf.org/html/rfc7946

The rupture format consists of a *FeatureCollection*, containing one
to many Features. The FeatureCollection should contain a dictionary
called *metadata*, which consists of the following fields:

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
| mech                  | (optional) Focal mechanism, one of "RS", "SS",        |
|                       | "NM", or "ALL".                                       |
+-----------------------+-------------------------------------------------------+

Each Feature should contain either a *Point* or *MultiPolygon*
geometry. When the geometry is a MultiPolygon, each sub-polygon in the
MultiPolygon must adhere to the following rules:

1. Each vertex should contain longitude, latitude *and* depth coordinates.
2. The vertices of the top edge are always defined first.
3. There must be the same number of vertices on top and bottom edges.
4. Last vertex must be the same as the first (polygon must be closed.)

While this is not required, if the top and bottom edges of each
sub-polygon are horizontal, the distance calculations are much faster.
   
The file should be named *rupture.json* and placed in the event's
*current* directory. A sample file is below::

  {
    "type": "FeatureCollection",
    "metadata": {
      "netid": "ci",
      "lon": -118.537,
      "lat": 34.213,
      "reference": " Source: Wald, D. J., T. H. Heaton, and K. W. Hudnut (1996). The Slip History of the 1994 Northridge, California, Earthquake Determined from Strong-Motion, Teleseismic, GPS, and Leveling Data, Bull. Seism. Soc. Am. 86, S49-S70.",
      "depth": 18.202,
      "time": "1994-01-17T12:30:55.390000Z",
      "locstring": "1km NNW of Reseda, CA",
      "mag": 6.7,
      "id": "ci3144585",
      "mech": "ALL",
      "rake": 0.0,
      "network": "California Integrated Seismic Network: Southern California Seismic Network (Caltech, USGS Pasadena, and Partners)"
    },
    "features": [
      {
        "type": "Feature",
        "properties": {
          "rupture type": "rupture extent"
        },
        "geometry": {
          "type": "MultiPolygon",
          "coordinates": [
          [
              [
                [
                  -118.421,
                  34.315,
                  5.0
                ],
                [
                  -118.587,
                  34.401,
                  5.0
                ],
                [
                  -118.693,
                  34.261,
                  20.427
                ],
                [
                  -118.527,
                  34.175,
                  20.427
                ],
                [
                  -118.421,
                  34.315,
                  5.0
                ]
              ]
            ]
          ]
        }
      }
    ]
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

