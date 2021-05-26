+---------+------------------+-----------------+------------+
| Azure   | |AzureM1015P37|  | OSX 10.15       | Python 3.7 |
+         +------------------+-----------------+------------+
|         | |AzureM1015P38|  | OSX 10.15       | Python 3.8 |
+         +------------------+-----------------+------------+
|         | |AzureM1014P37|  | OSX 10.14       | Python 3.7 |
+         +------------------+-----------------+------------+
|         | |AzureM1014P38|  | OSX 10.14       | Python 3.8 |
+         +------------------+-----------------+------------+
|         | |AzureLP37|      | ubuntu          | Python 3.7 |
+         +------------------+-----------------+------------+
|         | |AzureLP38|      | ubuntu          | Python 3.8 |
+---------+------------------+-----------------+------------+
| Travis  | |Travis|         | ubuntu          | Python 3.7 |
+---------+------------------+-----------------+------------+
| Codacy  | |Codacy|                                        |
+---------+-------------------------------------------------+
| CodeCov | |CodeCov|                                       |
+---------+-------------------------------------------------+


.. |Travis| image:: https://travis-ci.org/usgs/shakemap.svg?branch=master
    :target: https://travis-ci.org/usgs/shakemap
    :alt: Travis Build Status

.. |CodeCov| image:: https://codecov.io/gh/usgs/shakemap/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/usgs/shakemap
    :alt: Code Coverage Status

.. |Codacy| image:: https://api.codacy.com/project/badge/Grade/1f771008e85041b89b97b6d12d85298a
    :target: https://www.codacy.com/app/emthompson_2/shakemap?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=usgs/shakemap&amp;utm_campaign=Badge_Grade

.. |AzureM1015P37| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_10_15_Python37
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureM1015P38| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_10_15_Python38
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureM1014P37| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_10_14_Python37
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureM1014P38| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_10_14_Python38
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureLP37| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20Linux_Python37
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureLP38| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20Linux_Python38
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             



shakemap
========
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
