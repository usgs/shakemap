+---------------+----------------------+ 
| Linux build   | |Travis|             | 
+---------------+----------------------+ 
| Code quality  | |Codacy|             | 
+---------------+----------------------+ 
| Code coverage | |CodeCov|            | 
+---------------+----------------------+ 
| Azure builds  | |Azure|              | 
+---------------+----------------------+ 


.. |Travis| image:: https://travis-ci.org/vinceq-usgs/shakemap.svg?branch=altIntensity
    :target: https://travis-ci.org/vinceq-usgs/shakemap
    :alt: Travis Build Status

.. |CodeCov| image:: https://codecov.io/gh/vinceq-usgs/shakemap/branch/altIntensity/graph/badge.svg
    :target: https://codecov.io/gh/vinceq-usgs/shakemap
    :alt: Code Coverage Status

shakemap
========

branch: vinceq-usgs/shakemap/altIntensity

Note, this repository consists of code for version 4 of Shakemap.
The previous version (v3.5) of ShakeMap can be retrieved via 
Subversion (svn) from

- https://vault.gps.caltech.edu/repos/products/shakemap/tags/release-3.5/

The manual for ShakeMap V3.5 can be found here

- http://usgs.github.io/shakemap/manual3_5/index.html


Documentation
-------------

The ShakeMap v4 docs can be found `here <https://usgs.github.io/shakemap/sm4_index.html>`_.
Also, see the wiki associated with this repository for an install/setup
tutorial.


Introduction
------------

This repository holds the latest version of the ShakeMap code (v4.0).
ShakeMap is a system for rapidly characterizing the extent nad distribution of
strong ground shaking following significant earthquakes. The code is 
primarily written in the Python programming language. See the documentation 
at http://usgs.github.io/shakemap/sm4_index.html for a more detailed discussion
of ShakeMap and a list of references.

To receive updates on ShakeMap and discuss the software, please join the
`shake-dev <https://geohazards.usgs.gov/mailman/listinfo/shake-dev>`_
mailing list. We also encourage requests, questions, and discussions through
the Github issues tab associated with this repository.

Installation and Dependencies
-----------------------------

- Mac OSX or Linux operating systems
- bash shell, gcc, git, curl
- On OSX, Xcode and command line tools
- The ``install.sh`` script installs this package and all other dependencies,
  including python 3.5 and the required python libraries. It is regularly tested
  on OSX, CentOS, and Ubuntu.
- See the wiki for this repository for installation instructions.
