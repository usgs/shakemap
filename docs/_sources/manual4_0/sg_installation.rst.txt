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
software (ShakeMap runs on Python 3.7 or 3.8, while many systems still
use 2.7 as a default) and dependency headaches.

Configuration
=============

After running ``sm_profile`` the newly-created profile will have its 
*config* directory populated with a default set of configuration files.
These files should be reviewed by the operator prior to running 
ShakeMap. Additionally, the config files sometimes change when the
code is updated and so it is fairly common for errors encountered after
an update to be related to changes in the configs. We hope that the
configs will become more stable as the code matures, however it is
useful to check your configuration files against those in the
directory *shakemap/data* to make sure that you make any necessary 
changes to parameters.

Some internationalization support is included in the making of maps with
the ``mapping`` module. See the *products.conf* section below.

The configuration files are discussed in the sections below.

model.conf
----------

*model.conf* contains configuration options for the way modeling works, such
as the GMPE or GMPEs to use, the location of the Vs30, the IMTs to compute,
and options on where the predictions should be computed (i.e., grid
resolution or a list of site locations). One can make a copy of this
file in an event directory to have event-specific config options. 
In the event-specific *model.conf* it is only necessary to list parameters
that differ from those in the global file. Note that you must include
any section or sub-section indicators for the parameter in question. For
instance, to set the ``max_deviation`` parameter in an event-specific
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

*select.conf* has options for GMPE selection, which are used by
the ``select`` module. Note that if/when the ``select`` module runs, it
creates the file *model_select.conf* in the event's *current* directory,
which overrides the GMPE set in the *model.conf* file located in the
global config directory, but the config settings in an event-specific
*model.conf* take precedence over the settings in *model_select.conf*.
Thus, if there are any event-specific changes to the *model.conf*,
a sensible approach is to rename *model_select.conf* to *model.conf*
and then add any other config options to it.


Please see the
:ref:`Ground Motion Selection section <sec-select-4>` for
additional details on how this configuration works.


products.conf
-------------

*products.conf* holds options for the various ShakeMap products, such as
contours, rasters, and maps. Additional explanation is
available as comments in the *products.conf* file.

Of particular 
interest to some users will be the options for **language** and 
possibly **fontfamily**. These options allow for the maps to be
generated in a language other than English. If your language
of interest already exists in *shakemap/data/mapping* then you
can specify its extension for the **language** parameter. You may or
may not have to set the **fontfamily** parameter to support your
language of choice -- that will depend on your system configuration.

If your language does not appear as one of the options in 
*shakemap/data/mapping*, you may make the translations as
described in one of the language files (e.g., *map_strings.en*)
and add your new file (with the proper extension -- see
https://www.iana.org/assignments/language-subtag-registry
for a list of options) to *shakemap/data/mapping* directory, and then 
use that extension for your **language** parameter. Please also send us
your file so that we can include it in the repository. That way it will
be there the next time you update or install the ShakeMap software.

Note that the translations only apply to the products of the ``mapping``
module.


gmpe_sets.conf
--------------

This file defines the GMPE sets that are available to be set in
*model.conf*. These sets can be as simple as a single GMPE with a
weight of 1.0. The GMPE sets can be selected directly from *model.conf*,
or a the custom GMPE set created by the ``select`` module can be
selected.


modules.conf
------------

*modules.conf* ontrols what GMPEs are available for constructing GMPE sets.
Generally, this only needs to be edited if you wish to use a GMPE that is not
currently imported. The GMPEs are imported
from the `OpenQuake Engine <https://github.com/gem/oq-engine>`_
`hazardlib <https://github.com/gem/oq-engine/tree/master/openquake/hazardlib>`_
library.


shake.conf
----------

This configuration file is only for very general configuration options
relating
to the operation of ``shake``. It allows the operator to configure additional
repositories of ShakeMap modules ("plugins," if you will). It also allows
the user to set the modules for automatic runs, using the parameter
**autorun_modules**. The
general idea is that shake can be run specifying specific modules like this::

  shake <event id> module1 module2

But since there are many modules and ``shake`` is often invoked via
automated processes, it is convenient to configure a list of
***utorun_modules** that will be used when no module is specified
on the command line, i.e., like this::

  shake <event id>



logging.conf
------------

*logging.conf* contains options for logging. Most users will likely not
need to modify this file unless they wish to change the format of the
messages, date/time stamps, or other logging behavior.

transfer.conf
-------------

*transfer.conf* controls the transfer of ShakeMap products to remote
systems via the
``transfer`` module. See the documentation within the file itself for
explanation of the available options.

migrate.conf
------------

This file holds parameters that determine how ShakeMap 3.5 data directories
are migrated to ShakeMap v4-compatible directories via the program
``sm_migrate``. This file allows the user to choose which OpenQuake
GMPE should be used in place of the ShakeMap GMPE previously used
for each event.


Downloading and Configuring Vs30 and Topography
===============================================

We provide three files available by FTP at 
ftp://hazards.cr.usgs.gov/shakemap:

* *global_vs30.grd* -- The 30 arcsecond resolution Vs30 data set for the entire globe.
* *topo_30sec.grd* -- The 30 arcsecond resolution topography data for the entire globe.
* *topo_15sec.grd* -- The 15 arcsecond resolution topography data for the entire globe.

By 'entire globe' we mean from 56 degrees south to 84 degrees north latitude.

Note that ``sm_profile`` allows the user to download the 30-arcsecond topo
and Vs30 files as part of the creation of a profile. If ``sm_profile`` is
called with the ``-a`` option, these files will be downloaded automatically
and the profile will be configured to use them.

If you did not have ``sm_profile`` download the grids, you have a choice
of 15 or 30 second resolution topography. 15 second data shows
more detail at small scales, but causes ShakeMap to take significantly
longer to make the various output maps. The ShakeMap system at the National
Earthquake Information Center uses the 30 second data. If you plan to use
the 15 second data, modify the topo file name below to topo_15sec.grd. 

Note that these files are somewhat large: the 30-second topo is 238 Mb, the
30-second Vs30 is 582 Mb, and the 15-second topo is 745 Mb.

To download the files, do::

    > mkdir [home]/shakemap_data
    > mkdir [home]/shakemap_data/vs30
    > mkdir [home]/shakemap_data/topo
    > cd [home]/shakemap_data/vs30
    > curl ftp://hazards.cr.usgs.gov/shakemap/global_vs30.grd -o global_vs30.grd
    > cd [home]/shakemap_data/topo
    > curl ftp://hazards.cr.usgs.gov/shakemap/topo_30sec.grd -o topo_30sec.grd

By default, the system is configured to find the Vs30 and topography files in 
the locations described above. To set the paths to other locations or file
names::

    > cd [home]/shakemap_profiles/[profile]/install/config

Modify *model.conf* to change the line::

    vs30file = <DATA_DIR>/vs30/global_vs30.grd

to the location of your Vs30 data. Similarly, edit *products.conf* and
change the line::

    topography = <DATA_DIR>/topo/topo_30sec.grd

to the path to your topography file. Note that ShakeMap completes
the macro ``<INSTALL_DIR>`` for the profile in question, but you may set 
the paths to any absolute path on your system.
