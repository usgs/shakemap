|Travis| |CodeCov|

.. |Travis| image:: https://travis-ci.org/usgs/shakemap.svg?branch=master
    :target: https://travis-ci.org/usgs/shakemap
    :alt: Travis Build Status
.. |CodeCov| image:: https://codecov.io/gh/usgs/shakemap/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/usgs/shakemap
    :alt: Code Coverage Status
.. image:: https://api.codacy.com/project/badge/Grade/1f771008e85041b89b97b6d12d85298a
    :target: https://www.codacy.com/app/emthompson_2/shakemap?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=usgs/shakemap&amp;utm_campaign=Badge_Grade
.. image:: https://badge.waffle.io/usgs/shakemap.svg?columns=all
    :target: https://waffle.io/usgs/shakemap
    :alt: 'Waffle.io - Columns and their card count'

shakemap
========
Note, this repository consists of experimental code that is still under
development for the next major update to Shakemap. The current version
of ShakeMap can be found here

- https://vault.gps.caltech.edu/repos/products/shakemap/tags/release-3.5/


Documentation
-------------

The API docs can be found `here <https://usgs.github.io/shakemap/>`_.
Also, see the wiki associated with this repository for an install/setup
tutorial.


Introduction
------------

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

To receive updates on ShakeMap and discuss the software, please join the
`shake-dev <https://geohazards.usgs.gov/mailman/listinfo/shake-dev>`_
mailing list. We also encourage requests, questions, and discussions through
the Github issues tab associated with this repository.

Installation and Dependencies
------------

- Mac OSX or Linux operating systems
- Python 3
- Python libraries: numpy, scipy, rasterio, fiona, xlrd, pandas, basemap,
  shapely, h5py, gdal, descartes, openquake.engine, neicio,
  MapIO, matplotlib, pytables, lxml
- The ``install.sh`` script installs this package and all dependencies. It is
  regularly tested on OSX and Ubuntu.
- Note that for development, we have found that the ``-e`` for the pip install
  is useful.

OQ Hazard Library
-----------------

One of the significant factors driving the rewrite of ShakeMap into the Python
language was the availability of the library of Ground Motion Prediction
Equations (GMPEs) and other tools incorporated into the OpenQuake (OQ_)
Hazard Library (openquake.hazardlib_).
The OQ hazard library provided us with a broad range of
well-tested, high performance, open source global GMPEs. Due to constraints
imposed by the software architecture of earlier implementations of ShakeMap, the
development and validation of GMPE modules is time consuming and difficult, which
restricted the quantity and timeliness of the available modules. The OQ Hazard Library
provides a broad array of current GMPE and related hazard modules, as well as a
framework for easily adding new modules (whether by GEM or ShakeMap staff),
jumpstarting our efforts to re-implement ShakeMap.

.. _OQ: https://github.com/gem/oq-engine/#openquake-engine
.. _openquake.hazardlib: http://docs.openquake.org/oq-engine/stable/openquake.hazardlib.html
