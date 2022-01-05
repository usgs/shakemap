+---------+------------------+-----------------+-------------+
| Azure   | |AzureLP38|      | ubuntu          | Python 3.8  |
+         +------------------+-----------------+-------------+
|         | |AzureLP39|      | ubuntu          | Python 3.9  |
+         +------------------+-----------------+-------------+
|         | |AzureLP310|     | ubuntu          | Python 3.10 |
+         +------------------+-----------------+-------------+
|         | |AzureM11P38|    | OSX 11          | Python 3.8  |
+         +------------------+-----------------+-------------+
|         | |AzureM11P39|    | OSX 11          | Python 3.9  |
+         +------------------+-----------------+-------------+
|         | |AzureM11P310|   | OSX 11          | Python 3.10 |
+         +------------------+-----------------+-------------+
|         | |AzureM1015P38|  | OSX 10.15       | Python 3.8  |
+         +------------------+-----------------+-------------+
|         | |AzureM1015P39|  | OSX 10.15       | Python 3.9  |
+         +------------------+-----------------+-------------+
|         | |AzureM1015P310| | OSX 10.15       | Python 3.10 |
+---------+------------------+-----------------+-------------+
| CodeCov | |CodeCov|                                        |
+---------+--------------------------------------------------+


.. |CodeCov| image:: https://codecov.io/gh/usgs/shakemap/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/usgs/shakemap
    :alt: Code Coverage Status

.. |AzureLP38| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20Linux_Python38
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureLP39| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20Linux_Python39
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureLP310| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20Linux_Python310
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureM11P38| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_11_Python38
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureM11P39| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_11_Python39
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureM11P310| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_11_Python310
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureM1015P38| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_10_15_Python38
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureM1015P39| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_10_15_Python39
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             

.. |AzureM1015P310| image:: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_apis/build/status/usgs.shakemap?branchName=master&jobName=ShakeMap&configuration=ShakeMap%20MacOS_10_15_Python310
   :target: https://dev.azure.com/GHSC-ESI/USGS-ShakeMap/_build/latest?definitionId=2&branchName=master
   :alt: Azure DevOps Build Status                                             



.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black


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
ShakeMap is a system for rapidly characterizing the extent and distribution of
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
- bash shell, git, curl
- On OSX, Xcode and command line tools
- The ``install.sh`` script installs this package and all other dependencies,
  including python 3.X and the required python libraries. It is regularly tested
  on OSX, CentOS, and Ubuntu.
- See the wiki for this repository for installation instructions.
