shakemap
========

|Travis| |CodeCov|
|QuantifiedCode| |Codacy|

.. |Travis| image:: https://travis-ci.org/usgs/shakemap.svg?branch=master
    :target: https://travis-ci.org/usgs/shakemap
    :alt: Travis Build Status
.. |CodeCov| image:: https://codecov.io/gh/usgs/shakemap/branch/master/graph/badge.svg)
    :target:https://codecov.io/gh/usgs/shakemap
    :alt: Code Coverage Status
.. |QuantifiedCode| image::https://www.quantifiedcode.com/api/v1/project/b5c9b23475e443aa99806921f729f55a/badge.svg
    :target:https://www.quantifiedcode.com/app/project/b5c9b23475e443aa99806921f729f55a
    :alt: Code issues
.. |Codacy| image::https://api.codacy.com/project/badge/Grade/1d3b94ef3793456f861c67d7905e7de7
    :target:https://www.codacy.com/app/emthompson/shakemap?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=usgs/shakemap&amp;utm_campaign=Badge_Grade
    :alt: Codacy code issues

## Disclaimer
This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software. 

## Introduction

ShakeMap is a system for rapidly characterizing the extent nad distribution of
strong ground shaking following significant earthquakes. The current stable
production version of ShakeMap (V3.5) is largely written in Perl, but also
makes use of GMT (Generic Mapping Tools), MySQL, and many other programs.
The manual for ShakeMap V3.5 can be found here

- http://usgs.github.io/shakemap/

The manual includes an both a Technical Guide and and User's Guide. Installation
instructions can be found in the chapter titled Software & Implementation Guide. 

This repository hosts the ShakeMap V3.5 manual and is
also where we refactoring the code base into Python. The core ShakeMap code,
approaching fifteen years old, was overdue for a major overhaul to more
organically incorporate (or eliminate) the many extensions that had been added
over its lifetime, and to facilitate several new demands from ShakeMap’s
expanded role as a global provider of post-earthquake information and earthquake
scenarios, and as the input to loss modeling software.

ShakeMap was originally written for use at the Southern California Seismic
Network. Over time, it has been adopted by many national and international
seismic networks as the hazard mapping tool of choice. It is now in operation
at all regional seismic networks within the United States, and the Global
ShakeMap System at the USGS’s National Earthquake Information Center in Golden,
Colorado. The varied nature of its national and international installations has
required extensive modifications to the original source code. Additional uses of
ShakeMap, such as for scenario earthquakes and the ShakeMap Atlas, have also
required ongoing modification of the code. 

## Dependencies

- Mac OSX or Linux operating systems
- Python 3
- Python libraries: numpy, scipy, rasterio, fiona, xlrd, pandas, basemap,
  basemap-data-hires, shapely, h5py, gdal, descartes, oq-hazlib, neicio,
  MapIO, matplotlib, jupyter

## OQ Hazard Library

One of the significant factors driving the rewrite of ShakeMap into the Python
language was the availability of the library of Ground Motion Prediction
Equations (GMPEs) and other tools incorporated into the OpenQuake
([OQ](www.globalquakemodel.org/openquake/about/))
hazard library ([oq-hazardlib](github.com/gem/oq-hazardlib)).
The OQ hazard library provided us with a broad range of
well-tested, high performance, open source global GMPEs. Due to constraints
imposed by the software architecture of earlier implementations of ShakeMap, the
development and validation of GMPE modules is time consuming and difficult, which
restricted the quantity and timeliness of the available modules. The oq-hazardlib
provides a broad array of current GMPE and related hazard modules, as well as a
framework for easily adding new modules (whether by GEM or ShakeMap staff),
jumpstarting our efforts to re-implement ShakeMap.


