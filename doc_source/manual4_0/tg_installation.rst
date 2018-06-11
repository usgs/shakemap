.. _sec-installation-4:

******************************
Installation and Configuration
******************************

The `ShakeMap Wiki <https://github.com/usgs/shakemap/wiki>`_ provides
a basic quick-start guide to installing and running ShakeMap v4. The
present section is intended to provide supplementary material and
additional detail for installing, configuring, and running v4.

Installation
============

The Wiki does a pretty good job of explaining the installation process,
which is fairly automated. Here we will just reiterate that things will
go most smoothly if you use the bash shell and conda virtual environment.
Taking a more manual approach will likely lead to conflicts with system
software (ShakeMap runs on Python 3.6, while most systems still use 2.7
as a default) and dependency headaches.

Configuration
=============

After running **sm_profile** the newly-created profile will have its 
*config* directory populated with a default set of configuration files.
These files should be reviewed by the operator prior to running 
ShakeMap. The configuration files are discussed in the sections below:

model.conf
----------

Config options for how modeling works, such as what 
GMPE or GMPEs to use, what Vs30 file to use, what IMTs to compute, and
options on where the predictions should be computed (i.e., grid
resolution or a list of site locations). One can make a copy of this
file in an event directory to have event-specific config options. 
In the event-specific *model.conf* it is only necessary to list parameters
that differ from those in the global file. Note that you must include
any section or sub-section indicators for the parameter in question. For
instance, to set the `max_deviation` parameter in an event-specific
model.conf file, you would include the lines::

    [data]
      [[outlier]]
        max_deviation = 2.0

One may also create a custom GMPE set in the event-specific *model.conf*
and set the system to use it. For instance::

    [gmpe_sets]
      [[gmpe_northridge_custom]]
        gmpes = active_crustal_california,
        weights = 1.0,
        weights_larage_dist = None
        dist_cutoff = nan
        site_gmpes = None
        weights_site_gmpes = None
    [modeling]
      gmpe = gmpe_northridge_custom


select.conf
-----------

Config options for GMPE selection, which are used by
the `select` module. Note that if/when the `select` module runs, it
creates the file `model_zc.conf` in the event's _current_ directory,
which overrides the GMPE set in the `model.conf` file located in the
global config directory, but the config settings in an event-specific
`model.conf` take precedence over the settings in `model_zc.conf`.


products.conf
-------------

Options for the various ShakeMap products, such as
contours, rasters, and maps.


gmpe_sets.conf
--------------

This file defines the GMPE sets that are available to be set in
`model.conf`. These sets can be as simple as a single GMPE with a
weight of 1.0. The GMPE sets can be
selected directly in `model.conf`, or a new GMPE set is created by the
**select** module, which combines multiple GMPE sets associated with
specific tectonic environments to account for uncertainty in the
tectonic classification (e.g., if an event is located near the boundary
between active and stable regions).


modules.conf
------------

Controls what GMPEs are available for
constructing GMPE sets. Generally, this only needs to be edited if you
wish to use a GMPE that is not currently imported. The GMPEs are imported
from the `OpenQuake Engine <https://github.com/gem/oq-engine>`_
`hazardlib <https://github.com/gem/oq-engine/tree/master/openquake/hazardlib>`_
library.


shake.conf
----------

This config file is only for very general configuration options relating
to the operation of **shake**.


logging.conf
------------

Contains options for logging. Most users will likely not need to modify
this file unless they wish to change the format of the messages or other
logging behavior.

transfer.conf
-------------

Controls the transfer of ShakeMap products to remote systems via the
**transfer** module. See the 
documentation within the file itself for explanation of the available
options.

migrate.conf
------------

Parameters that determine how ShakeMap 3.5 data directories are 
migrated to ShakeMap 4.0-compatible directories via the program
**sm_migrate**. This file allows the user to choose which OpenQuake
GMPE should be used in place of the ShakeMap GMPE previously used
for each event.


Downloading and Configuring Vs30 and Topography
===============================================

We provide three files available by FTP at 
ftp://hazards.cr.usgs.gov/shakemap:

* `global_vs30.grd` -- The 30 arcsecond resolution Vs30 data set for the entire globe.
* `topo_30sec.grd` -- The 30 arcsecond resolution topography data for the entire globe.
* `topo_15sec.grd` -- The 15 arcsecond resolution topography data for the entire globe.

By 'entire globe' we mean from 56 degrees south to 84 degrees north latitude.

You have a choice of 15 or 30 second resolution topography. 15 second data shows
more detail at small scales, but causes ShakeMap to take *significantly* longer
to make the various output maps. The ShakeMap system at the National Earthquake
Information Center uses the 30 second data. If you plan to use the 15 second
data, modify the topo file name below to topo_15sec.grd. 

To download the files, do::

    > cd [home]/ShakeMap/[profile]/install/data/mapping

    > curl ftp://hazards.cr.usgs.gov/shakemap/global_vs30.grd -o global_vs30.grd

    > curl ftp://hazards.cr.usgs.gov/shakemap/topo_30sec.grd -o topo_30sec.grd

By default, the system is configured to find the Vs30 and topography files in 
the locations described above. To set the paths to other locations or file
names::

    > cd [home]/ShakeMap/[profile]/install/config

Modify `model.conf` to change the line::

    vs30file = <INSTALL_DIR>/data/mapping/global_vs30.grd

to the location of your Vs30 data. Similarly, modify products.conf to
change the line::

    topography = <INSTALL_DIR>/data/mapping/topo_30sec.grd

to the path to your topography file. Note that ShakeMap completes
the macro <INSTALL_DIR> for the profile in question, but you may set 
the paths to any absolute path on your system.
