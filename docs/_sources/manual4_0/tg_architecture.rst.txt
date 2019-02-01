.. _sec-architecture-4:

****************************
Software Architecture
****************************

Overview
========

ShakeMap 4.0 is designed to allow flexibility in the organization of
computing resources. :num:`Figure #architecture-overview` illustrates an
idealized implementation where data preparation, processing, and rendering
all take place within separate computational units. The processing sequence
starts when an earthquake is identified and a decision is made to produce a
ShakeMap. The process ``shake assemble`` collects the available information
about the event (origin and rupture parameters, seismic data, etc.) as well
as ShakeMap configuration information (which may include information about the
event's seismotectonic regime and related choices about GMPE selection), and
produces a file, *shake_data.hdf*, containing all of these parameters. This
file may be injected into a messaging system, but may also be used locally
by subsequent processes.

.. _architecture-overview:

.. figure:: _static/sm4.*
   :width: 700
   :align: left

   Overview of ShakeMap architecture.

The processing continues when *shake_data.hdf* becomes available. The ShakeMap
process ``shake model`` reads *shake_data.hdf* and produces output in the file
*shake_result.hdf*. This result can then be fed into a messaging system for
delivery to consumers. Some consumers, however, have more sophisticated
requirements than can be accommodated by simply processing *shake_result.hdf*
or other generic ShakeMap products.
ShakeCast :ref:`(Wald et al., 2008) <wald2008shakecast>`, for example, requires
ground motions at a variety of spectral periods and at specific locations that
may not fall on or within the grid produced by the authoritative ShakeMap
system. ShakeCast operators may also have data not available to the
authoritative system. Remote processing systems can receive *shake_data.hdf*
from a messaging system, and run the program ``shake augment`` to add their own
data and configuration choices to those contained in *shake_data.hdf*
(see :num:`Figure #shake-consumer`). They may then run ``shake model`` to
generate a *shake_result.hdf* specific to their needs.

.. _shake-consumer:

.. figure:: _static/consumer.*
   :width: 700px
   :align: left

   An example of a consumer of the *shake_data.hdf* product.

Rendering begins when *shake_result.hdf* becomes available. A set of modules
exist (or may be developed) to read *shake_result.hdf* and produce the variety
of products for which ShakeMap is known.
These ``shake`` modules may produce the products locally (i.e., by the same
system that generates *shake_result.hdf*) and transfer them to consumers via
a messaging system or other means.

An alternative approach, however, is to create a web service that delivers
the products when they are requested. This approach is illustrated in
:num:`Figure #shake-web`. When the website is notified of the existence
of *shake_result.hdf*, it can begin the process of creating a "page" for the
event. It requests any necessary products from the web service, which in turn
generates those products from *shake_result.hdf* (via ``shake`` modules). As
products are needed (e.g., from users viewing or requesting downloads) they
are produced on the fly by the web service. Once generated, products may be
cached by the web system to improve performance.

.. _shake-web:

.. figure:: _static/web.*
   :width: 700px
   :align: left

   An example of a website using a web service to retrieve products. The web
   service produces products from *shake_result.hdf*.

Any combination of these approaches (i.e., producing products locally or
via a web service) may be developed (e.g., the web service may be designed
to collect a subset of ShakeMap products available through a messaging
system and deliver them when requested, rather than producing them itself).
Thus, the same set of constituent modules are needed, whether the products
are delivered directly by the authoritative ShakeMap system or through a
web service.

Programs
========

The core components of ShakeMap are a set of command line programs.
These programs allow the operator to set up a ShakeMap environment,
collect data and configurations into inputs (i.e., *shake_data.hdf*),
and generate ShakeMap grids and their associated products.

sm_profile
----------

The user will need to run ``sm_profile`` at least once to create a ShakeMap
environment, referred to as a 'profile.' This environment consists of two
directories -- one for event data, and another for configuration files and
associated support products (Vs30 grid, geographic features, etc.) -- and a
configuration file that points to them. The profile data resides in a file
called *profiles.conf* in a subdirectory, *.shakemap*, of the user's home
directory. The user may choose another location for the profile file by using
the ``-f`` option to ``sm_profile``. Other ShakeMap programs read the profile
information and use it to find event and configuration information. By default,
the install and data paths are located in a directory
``[home]/shake_profiles/[profile]`` where
``[home]`` is the user's home directory, and ``[profile]`` is the name of
the newly-created profile (often "default"). It is advisable, however
to make the data path point to a disk with a large amount of free space
(via either directly setting it when creating the profile, or through the
use of symbolic links).

The data directory (``<data_dir>``) contains event subdirectories (named
with their event IDs) and their associated subdirectories::

    <data_dir>/
        <event_id_1>/
            current/
                event.xml
                *_dat.xml
                *_fault.txt (or rupture.json)
                model.conf (or model_select.conf)
                products/
                    shake_result.hdf
                    ...
            shake.log
            .backup0001/
                event.xml
                ...
            .backup0002/
                ...
            ...
        <event_id_2>/
            ...
        <event_id_3>/
            ...
        ...

The 'install' directory (``<install_dir>``) holds configuration files and
user supplied geographic or other system specific data::

    <install_dir>/
        config/
            model.conf
            modules.conf
            gmpe_sets.conf
            ...
        mapping/
            global_vs30.grd
            global_topo.grd
        logs/
            shake.log
            shake.log.<timestamp>
            ...
        <other_directory>/
            (additional data files)
        ...

Macros within the configuration system allow the user to specify the
root data and install directories when setting configuration
parameters.

The user may have more than one profile, and can switch between them with
``sm_profile``. This allows the user to have different configurations
and data repositories for different event sets (e.g., real time events,
scenarios, and historic events). See the
:ref:`sm_profile man page <sm-profile>` for usage and a list of options.

shake
-------

The primary ShakeMap program is called ``shake``. It takes an event ID and an
optional list of modules as arguments. The modules do the work of assembling
the input data, producing interpolated grids, and deriving products from the
grids and associated metadata. See the :ref:`shake man page <shake>` or run
``shake --help`` for a list of available modules. Each of the modules may have
its own command-line options; run ``shake help MODULE`` to see the help for a
given module. 

The behavior of ``shake`` and some of its modules are controlled by
the configuration files *shake.conf*, *logging.conf* and *products.conf*.
*logging.conf* is largely concerned with the way logging is handled, it has 
relatively few parameters that the users will want to change. 
*shake.conf* allows the user to specify alternative locations to
search for ShakeMap modules, and to configure a list of modules that will
be run when ``shake`` is called without any modules (this is useful as a
stanard set of modules that are run routinely). 
*products.conf* controls the behavior of some of the core modules that
produce ShakeMap products. See the comments in the config files for
additional details.
All three files should be in the user's current profile's
``INSTALL_DIR/config`` directory.

``shake`` will attempt to assess if a module is being run out of sequence
or if its dependencies are out of date. For instance, if a new data file
has been added, and the user tries to run the ``model`` module before
re-running the ``assemble`` module, ``shake`` will warn the user and 
quit without running ``model``. The user can override this behavior by
calling ``shake`` with the ``--force`` option.

shake Modules
-------------

Below is a description of many of the modules available to ``shake``.
They are ordered in more or less the order they would be called. The
module ``select`` would be run first if the operator wanted to have
the ShakeMap system determine the GMPE set, IPE, GMICE, or
cross-correlation functions to use based on the event's location and depth.
Many operators will have a fixed configuration for their regions, and will
therefore not use ``select``.
The operator will then typically run ``assemble`` (or possibly ``augment``)
to
create (or update) the *shake_data.hdf* input file, followed by running
``model``. The order of modules after ``model`` is usually not as important
as they typically depend upon the output of ``model`` (i.e.,
*shake_result.hdf*) and not upon the results of other modules.
``transfer``, however, will typically be run last.

select
```````
Please see the
:ref:`Ground Motion Selection section <sec-select-4>` for
details on how ``select`` works.


dyfi
````

The ``dyfi`` module queries ComCat for any "Did You Feel It?" data associated
with an event and writes that data to a file in the event's *current*
directory.
The event ID must be an ID that the ComCat system recognizes, thus the use of
event IDs other than those produced by NEIC or the US regional networks is
unlikely to produce results.

See :meth:`shakemap.coremods.dyfi` for the module's API documentation.

assemble
````````

The ``assemble`` module collects event and configuration data and creates the
file *shake_data.hdf*. It first reads *event.xml* and stores it in a data
structure. ``sm_assemble`` then reads the configuration files

.. code-block:: python

    <install_dir>/modules.conf
    <install_dir>/gmpe_sets.conf
    <install_dir>/model.conf


and assembles them into a single configuration. It then reads and 
incorporates

.. code-block:: python

    <data_dir>/<evnt_id>/current>/model.conf (or model_select.conf).

Any parameter set in the event-specific *model.conf* will override
parameters set in the other configuration files. Note: if both
*model.conf* and *model_select.conf* exist in the event directory,
*model.conf* will be processed and *model_select.conf* will be ignored.

``assemble`` then reads any files with a *_dat.xml* extension
and assembles them into a station list. See ??? for a description
of the data file format. Similarly, ``assmeble`` will read a
file with the *_fault.txt* (or *_fault.json*) extension and
process it as a specification of a finite rupture. See ??? for
a description of the rupture file formats. Note: only one rupture
file should be present in the event's input directory. If more
than one file exists, only the first (lexicographically) will we
processed.

If no backups exist (i.e., event subdirectories named *.backup????*)
then the ShakeMap history from an existing *shake_data.hdf* is
extracted and updated. If there is no current *shake_data.hdf*, the
history for the event is initiated. If backups do exist, then the
history is extracted from the most current backup and appended
with the current timestamp, originator, and version.

``assemble`` then consolidated all of this data and writes
*shake_data.hdf* in the event's *current* directory. If *shake_data.hdf*
already exists in that location, it will be overwritten.

``assemble`` takes an optional command-line argument (``-c COMMENT``
or ``--comment COMMENT``) to provide a comment
that will be added to the history for the
current version of the event's ShakeMap. If run from a terminal,
and a comment is not provided on the command line, ``assemble`` 
will prompt the user for a comment.
Run ``shake help assemble`` for more.

See :meth:`shakemap.coremods.assemble` for the module's API
documentation.

.. _shake-assemble:

.. figure:: _static/assemble.*
   :width: 700px
   :align: left

   Data flow of the *assemble* module.

augment
```````

The ``augment`` module behaves in a manner similar to ``assemble` except
that it will first read *shake_data.hdf* from the event's *current*
directory. If *exven.xml* exists in the event's *current* directory, its
data will replace the data in the existing *shake_data.hdf*.

The configuration data in *shake_data.hdf* is used as a starting point,
and any configuration data from the system configuration files or the
event's *model.conf* (or *model_select.conf*) will then be added to it. Where
there are conflicts, the system configuration parameters will override
those found in *shake_data.hdf*. The event-specific configuration
parameters from the local system retain the highest priority.

Data files (i.e., files in the event's *current* directory that have
the *_dat.xml* extension) will be added to any data already found in
*shake_data.hdf*. If a fault file is found in the local directory, it
will replace the existing fault data in *shake_data.hdf*.

The history information will be updated to reflect the update time and
originator (if applicable).

As with ``assemble``, ``augment`` takes an optional command-line 
argument (``-c COMMENT``
or ``--comment COMMENT``) to provide a comment
that will be added to the history for the
current version of the event's ShakeMap. If run from a terminal,
and a comment is not provided on the command line, ``assemble`` 
will prompt the user for a comment.
Run ``shake help augment`` for more.

See :meth:`shakemap.coremods.augment` for the module's API 
documentation.

model
`````

The ``model`` module reads the data in *shake_data.hdf* and produces an
interpolated ShakeMap. Depending upon the settings found in *model.conf*,
the interpolation product may be a grid or a set of points. See
*model.conf* for additional options and documentation. The *model.conf*
file in the user's current profile (i.e., ``INSTALL_DIR/config/model.conf``)
will be read first, and then if *model.conf* or *model_select.conf* exists
in the event's *current* directory, then the parameters set therein will
override those in the profile's *model.conf*. If both *model.conf* and
*model_select.conf* exist in the event's *current* directory, *model.conf*
will be read and *model_select.conf* will be ignored. ``model`` also reads
the configuration files *gmpe_sets.conf* and *modules.conf*, which
reside in the current profile's ``INSTALL_DIR/config`` directory. See
the documentation within those files for more information.

A great deal of this manual is devoted to the way the interpolation is
performed, and the effect of various configuration options. See the
relevant sections for more. In particular, the section
:ref:`sec-processing-4`
goes into detail on the way the model program works.

``model`` writes a file, *shake_result.hdf*, in the *products*
subdirectory of the event's *current* directory.
See :ref:`the formats section <sec-output-formats-4>`
of this manual for more on the format and
content of *shake_result.hdf*.

See :meth:`shakemap.coremods.model` for the module's API
documentation.

contour
```````

``contour`` reads an event's *shake_result.hdf* and produces iso-seismal
contours for each of the intensity measure types found therein. The contours
are written as GeoJSON to files called *<imt_type>_cont.json* in the event's
*current/products* subdirectory.

See :meth:`shakemap.coremods.contour` for the module's API documentation..

gridxml
```````

``gridxml`` reads an event's *shake_result.hdf* and produces the ShakeMap 3.5
files *grid.xml* and *uncertainty.xml*. Note that these files will eventually
become deprecated in favor of the new *shake_result.hdf* file.
See the products section of this manual for more on these files.
Note that the use of these files is
deprecated. System designers should extract the relevant information
directly from *shake_result.hdf*. See :ref:`the formats section <sec-output-formats-4>`
of this manual for more on
using *shake_result.hdf*.

See :meth:`shakemap.coremods.gridxml` for the module's API
documentation.

info
```````

``info`` reads an event's *shake_result.hdf* and produces *info.json*,
which contains metadata about the ShakeMap.

See :meth:`shakemap.coremods.info` for the module's API
documentation.

kml
````

``kml`` reads *shake_result.hdf* and produces a file, *shakemap.kmz*
which is a self-contained file of geographic layers suitable for 
reading into Google Earth. The layers include an MMI overlay, MMI
contours, station locations and data, and the event's epicenter.

mapping
```````

``mapping`` reads an event's *shake_result.hdf* and produces a set of
maps of the IMTs.

See :meth:`shakemap.coremods.mapping` for the module's API documentation.
See the configuration file *products.conf* for information on configuring
the ``mapping`` module.

plotregr
````````

``plotregr`` makes plots of an approximation of the GMPE's predicted
ground motion as a function of distance for each output IMT, along with
the data for that IMT. The plotted value at each distance is the average
value of the GMPE's gridded values in that particular distance bin. The
values are predicted on soil (Vs30=180 m/s) and rock (Vs30=760 m/s), which
are plotted as green and red lines, respectively. The +/- 1 standard
deviation lines are also plotted. The station and dyfi data are plotted at
their computed distance from the source. If no finite fault is available for
the map, then the approximated point-source to finite-fault  distance 
is used.

The ``plotregr`` module is fairly simplistic and of limited utility. 
Our hope it that it will be rendered obsolete by
more sophisticated interactive plots.

raster
```````

``raster`` reads an event's *shake_result.hdf* and produces GIS
raster files of the mean and standard deviation for each of the
IMTs in *shake_result.hdf*.

See :meth:`shakemap.coremods.raster` for the module's API
documentation.

rupture
```````

``rupture`` reads an event's *shake_result.hdf* and produces a
file, *rupture.json* containing the coordinates of the rupture
plane(s) supplied via the input file *<>_fault.txt* or *<>_fault.json*.

See :meth:`shakemap.coremods.rupture` for the module's API
documentation.

shape
`````

``shape`` reads an event's *shake_result.hdf* and produces a set of
ESRI-style shape files for the ground motion parameters found therein.
The shape files are zipped together with supporting .lyr, .prj, and
metadata XML files and distributed as a file called *shape.zip*. The
use of the shape files is deprecated and it is preferable to use the
ESRI raster file (see the ``raster`` section, above.)

stations
````````

``stations`` reads an event's *shake_result.hdf* and produces a
JSON file, *stationlist.json*, of the input station data.

See :meth:`shakemap.coremods.stations` for the module's API
documentation.

transfer
`````````

``transfer`` allows the operator to transfer ShakeMap products to
other system via PDL or ssh. See the documentation in *transfer.conf*
for details on configuring ``transfer``.

See :meth:`shakemap.coremods.transfer` for the module's API
documentation.


Additional Programs
===================

ShakeMap provides a few auxiliary programs that may occasionally be 
useful.

getdyfi
-------

``getdyfi`` is a standalone program implementing the ``dyfi`` module's
functionality.  See the :ref:`getdyfi man page <getdyfi>` for usage and a 
list of options.

sm_compare
-----------

Allows the user to compare two ShakeMaps by making images of their
difference and ratio. 

See the :ref:`sm_compare man page <sm-compare>` for usage and a 
list of options.


sm_create
---------

``sm_create`` queries the NEIC ComCat database for ShakeMap data 
associated with an event and writes that data into the event's
local *current* directory. The event will then be available for
local processing.

See the :ref:`sm_create man page <sm-create>` for usage and a 
list of options.


sm_migrate
-----------

Migrates a directory of ShakeMap 3.5 data files into ShakeMap 4.0 inputs.
The migration of GMPEs is configurable via the *migrate.conf* configuration
file.

See the :ref:`sm_migrate man page <sm-migrate>` for usage and a 
list of options.


sm_queue
-----------

A daemon process to receive messages from external systems for the
triggering and cancellation of ShakeMap runs. See the section
:ref:`Queueing Events <sec-queue-4>` for more detail. The behavior
of ``sm_queue`` is controlled by the *queue.conf* configuration
file. This file is not copied to new profiles by default, so it 
may be retrieved from the source directory 
``<shake_install_dir>/shakemap/data``.

See the :ref:`sm_queue man page <sm-queue>` for usage and a 
list of options.


receive_amps, receive_origins, and associate_amps
-------------------------------------------------

``receive_amps`` and ``receive_origins`` are intended to be run 
by a configured instance of the USGS's Product
Distribution system to inform ``sm_queue`` of new origins and
unassociated amplitudes. They are, therefore, of limited utility
to most users, however they may serve as guides as to writing
similar programs for other systems.

``associate_amps`` is a utility program to associate the 
unassociated amplitudes with origins, and to create ShakeMap
input files with those that associate. Again, this will be of
limited utility to users not running PDL.


run_verification
----------------

Runs a set of simple verification tests and displays the results.
The resulting plots may be compared to those found in the 
documentation section :ref:`Verification <sec-verification-4>`.
``run_verification`` is a shell script. See the source file for
usage and notes.



