.. _sec-architecute-4:

****************************
Software Architecture
****************************

Overview
========

ShakeMap 4.0 is designed to allow flexibility in the organization of 
computing resources. :num:`Figure #architecture-overview` illustrates an 
idealized implementation where data preparation, processing, and rendering 
all take place within separate computational units. The processing 
sequence starts when an earthquake is identified and a decision is made
to produce a ShakeMap. The process **sm_assemble**
collects the available information about the event (origin and rupture 
parameters, seismic data, etc.) as well as ShakeMap configuration 
information (which may include information about the event's 
seismotectonic regime and related choices about GMPE selection), and 
produces a file, *shake_data.hdf*, containing all of these parameters. This 
file may be injected into a messaging system, but may also be used locally 
by subsequent processes. 

.. _architecture-overview:

.. figure:: _static/sm4.*
   :width: 650px
   :align: left

   Overview of ShakeMap architecture.

The processing continues when *shake_data.hdf* becomes available. The ShakeMap 
program **sm_model** reads *shake_data.hdf* and produces output in the file 
*shake_result.hdf*. This result can then be fed into a messaging system for 
delivery to consumers. Some consumers, however, have more sophisticated 
requirements than can be accommodated by simply processing *shake_result.hdf* 
or other generic ShakeMap products. 
ShakeCast :ref:`(Wald et al., 2008) <wald2008shakecast>`, for example, requires 
ground motions at a variety of spectral periods and at specific locations that 
may not fall on or within the grid produced by the authoritative ShakeMap 
system. ShakeCast operators may also have data not available to the 
authoritative system. Remote processing systems can receive *shake_data.hdf* 
from a messaging system, and run the program **sm_augment** to add their own 
data and configuration choices to those contained in *shake_data.hdf* 
(see :num:`Figure #shake-consumer`). They may then run **sm_model** to 
generate a *shake_result.hdf* specific to their needs.

.. _shake-consumer:

.. figure:: _static/consumer.*
   :width: 650px
   :align: left

   An example of a consumer of the *shake_data.hdf* product.

Rendering begins when *shake_result.hdf* becomes available. A set of programs 
can be developed to read *shake_result.hdf* and produce the variety of products 
for which ShakeMap is known. These programs may be wrapped together under the 
general title **sm_genex** (ShakeMap GENerate EXport). **sm_genex** may produce 
the products locally (i.e., by the same system that generates *shake_result.hdf*) 
and transfer them to consumers via a messaging system or other means.

An alternative approach, however, is to create a web service that delivers 
the products when they are requested. This approach is illustrated in 
:num:`Figure #shake-web`. When the website is notified of the existence 
of *shake_result.hdf*, it can begin the process of creating a "page" for the 
event. It requests any necessary products from the web service, which in turn 
generates those products from *shake_result.hdf* (via **sm_genex**). As 
products are needed (e.g., from users viewing or requesting downloads) they 
are produced on the fly by the web service. Once generated, products may be 
cached by the web system to improve performance.

.. _shake-web:

.. figure:: _static/web.*
   :width: 650px
   :align: left

   An example of a website using a web service to retrieve products. The web 
   service produces products from *shake_result.hdf*.

Any combination of these approaches (i.e., producing products locally or via a 
web service) may be developed (e.g., the web service may be designed to collect 
a subset of ShakeMap products available through a messaging system and deliver 
them when requested, rather than producing them itself). Thus, the same set of 
constituent modules are needed, whether the products are delivered directly by 
the authoritative ShakeMap system or through a web service.

Programs
========

The core components of ShakeMap are a set of command line programs. These 
programs allow the operator to set up a ShakeMap environment, collect 
data and configurations into inputs (i.e., *shake_data.hdf*), and
generate ShakeMap grids and their associated products. 

sm_profile
----------

The user will need to run **sm_assemble** at least once to create a 
ShakeMap
environment, referred to as a 'profile.' This environment consists of two 
directories -- one for 
event data, and another for configuration files and associated support
products (Vs30 grid, geographic information, etc.) -- and a configuration
file that points to them. The profile data resides in a file called
*profiles.conf* in a subdirectory, *.shakemap*, of the user's home 
directory. Other ShakeMap programs read the profile 
information in this configuration file and use it to find event and 
configuration information.

The data directory ('<data_dir>') contains event subdirectories (named
with their event IDs) and their associated subdirectories::

    <data_dir>/
        <event_id_1>/
            current/
                event.xml
                *_dat.xml
                *_fault.txt (or .json)
                model.conf (or model_zc.conf)
                products/
                    shake_result.hdf
                    ...
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

The 'install' directory ('<install_dir>') holds configuration files and 
user supplied geographic or other system specific data::

    <install_dir>/
        config/
            model.conf
            modules.conf
            gmpe_sets.conf
            ...
        site_data/
            vs30.grid
        <other_directory>/
            (additional data files)
        ...

Macros within the configuration system allow the user to specify the 
root data and install directories when setting configuration 
parameters.

The user may have more than one profile, and switch between them with
**sm_profile**. This allows the user to have different configurations 
and data repositories for different event sets (e.g., real time events,
scenarios, and historic events). See the 
:ref:`sm_profile man page <sm-profile>` for usage and a list of options.

sm_assemble
-----------

**sm_assemble** collects event and configuration data and creates the 
file *shake_data.hdf*. It first reads *event.xml* and stores it in a 
data structure. **sm_assemble** then reads the configuration files

| <install_dir>/modules.conf
| <install_dir>/gmpe_sets.conf
| <install_dir>/model.conf

and assembles them into a single configuration. It then reads 

| <data_dir>/<evnt_id>/current>/model.conf (or model_zc.conf).

Any parameter set in the event-specific *model.conf* will override 
parameters set in the other configuration files. Note: if both 
*model.conf* and *model_zc.conf* exist in the event directory, 
*model.conf* will be processed and *model_zc.conf* will be ignored.

**sm_assemble** then reads any files with a *_dat.xml* extension 
and assembles them into a station list. See ??? for a description 
of the data file format. Similarly, **sm_assmeble** will read a
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

**sm_assemble** then consolidated all of this data and writes 
*shake_data.hdf* in the event's *current* directory. If *shake_data.hdf*
already exists in that location, it will be overwritten.

See the :ref:`sm_assemble man page <sm-assemble>` for usage.

.. _shake-assemble:

.. figure:: _static/assemble.*
   :width: 650px
   :align: left

   Data flow of the *sm_assemble* process.

sm_augment
----------

**sm_augment** behaves very similar to **sm_assemble** except that it
will first read *shake_data.hdf* from the event's *current* directory.
If *exven.xml* exists in the event's *current* directory, its data will
replace the data in the existing *shake_data.hdf*. 

The configuration data in *shake_data.hdf* is used as a starting point,
and any configuration data from the system configuration files or the 
event's *model.conf* (or *model_zc.conf*) will then be added to it. Where
there are conflicts, the system configuration parameters will override 
those found in *shake_data.hdf*. The event-specific configuration 
parameters from the local system retain the highest priority.

Data files (i.e., files in the event's *current* directory that have 
the *_dat.xml* extension) will be added to any data already found in
*shake_data.hdf*. If a fault file is found in the local directory, it 
will replace the existing fault data in *shake_data.hdf*.

The history information will be updated to reflect the update time and
originator (if applicable).

See the :ref:`sm_augment man page <sm-augment>` for usage.

sm_model
--------

**sm_model** reads the data in *shake_data.hdf* and produces an
interpolated ShakeMap. Depending upon the settings found in *model.conf*,
the interpolation product may be a grid or a set of points. See
*model.conf* for additional options and documentation.

A great deal of this manual is devoted to the way the interpolation is
performed, and the effect of various configuration options. See the 
relevant sections for more.

**sm_model** writes a file, *shake_result.hdf*, in the *products*
subdirectory of the event's *current* directory. Aside from interpolated
grids (or lists of points) for each intensity measure type, and their
associated uncertainties, **sm_model** also writes ShakeMap metadata
('*info.json*'), station data ('*stationlist.json*'), the finite
fault information ('*rupture.json*'), and the entire configuration
object into *shake_result.hdf*. See ??? for more on the format and
content of *shake_result.hdf*.

See the :ref:`sm_model man page <sm-model>` for usage.

sm_contour
----------

**sm_contour** reads an event's *shake_result.hdf* and produces iso-seismal
contours for each of the intensity measure types found therein. The contours
are written as GeoJSON to files called *<imt_type>_cont.json* in the event's 
*current/products* subdirectory.

See the :ref:`sm_contour man page <sm-contour>` for usage.

sm_select
----------

**sm_select** reads an event's *event.xml* file for origin information
and then constructs a GMPE set for the event based on the event's residence 
within,
and proximity to, a set of predefined tectonic regions. The GMPE set, and the
selection of that GMPE set for use in processing, are written to 
*model_zc.conf* in the event's *current* directory.

The behavior of **sm_select** is controlled by the *select.conf*
configuration file. See the documentation in *select.conf* for more on
customizing **sm_select**.

The tectonic regions, and additional geographic layers, that the event
may fall within are defined by the STREC configuration. See the STREC
documentation for information on adding additional layers, then use
*select.conf* to customize the GMPE sets that the new layers will use.

The process by which sm_select builds a GMPE set is somewhat complicated.
STREC reports the tectonic region the earthquake lies within, as well
as the distance to the closest polygon of the other tectonic region
types. For example, for an earthquake in California STREC would report
that the event was zero distance from region 'acr'
(which is to say that it lies within the active crustal region), but
STREC would also report distances to regions 'scr' (stable continental),
'volcanic', and 'subduction'. Each non-subduction region is also 
configured with a "horizontal buffer." The buffer determines how far
the region extends into neighboring regions. The buffer for subduction
regions is always zero.\ [#fn1]_ If the event happens within the buffer
of a neighboring region, the distance and buffer are used to build a 
weighted combination of the GMPE sets representing the regions in 
question.

For example, if an earthquake occurred within the 'scr' region, but
was 40 km from the "acr" region, and the 'acr' region's horizontal
buffer was 100 km, then the 'scr' region would be given a weight
of 1.0, and the 'acr' would be given (100 - 40) / 100 = 0.6. Normalizing
by the total, the final weights would be 0.625 'scr' and 0.375 'acr'.

Each region's GMPE set is in turn comprised of a weighted set of other
GMPE sets, based on the earthquake's depth. For each of the non-subduction
regions,
*sm_select* builds a weighted combination of the configured GMPE sets
based on the event's depth. If the earthquake falls within a subduction 
regions, STREC 
reports the probabilities that the earthquake is crustal, on the 
subduction interface, or within the subducting slab. *sm_select* 
combines the GMPE sets for each of these regimes, weighted by their
probabilities, into a subduction GMPE set that is specific to the 
earthquake's location.

*sm_select* also considers the earthquake's presence within, or 
distance from, 
any number of user-defined geographic layers. If the earthquake is 
within a layer, that layer's 
parameters (as configured in *select.conf*) replace the any or all
of the parameters of the corresponding tectonic regions, and the
calculation of a weighted GMPE set proceeds as before. For example,
the layer section of *select.conf* might contain::

    [layers]
        [[california]]
            horizontal_buffer = 50
            [[[scr]]
                horizontal_buffer = 25
            [[[acr]]]
                horizontal_buffer = 25
                gmpe = Special_California_GMPE
                min_depth = -Inf
                max_depth = Inf

If an earthquake falls within the 'california' layer, the tectonc regions 
'scr' and 'acr' would have their horizontal buffers reset to 25 km and,
in addition, the 'acr' region would have its GMPE selection reset to the 
GMPE set 'Special_California_GMPE' for earthquakes of all depths.

If the 
earthquake is not inside a custom geographic layer, but within the horizontal
buffer distance of one, the
GMPE sets for the modified and unmodified tectonic regions are each
determined separately
and a weighted combination of the two is computed (where the weights 
are based on the distance and the horizontal buffer, as described
above).

Unlike the tectonic regions, with the geographic layers only the nearest 
layer is considered. If an earthquake falls
within more than one layer (possible if layers are nested), the first one
encountered in the *select.conf* is used and any other(s) will be ignored.

See the :ref:`sm_select man page <sm-select>` for usage.


.. rubric:: Footnotes

.. [#fn1] Subduction regions do not extend beyond their defined boundaries
   because STREC cannot provide the crustal, interface, 
   and slab probabilities for events outside of the defined subduction 
   zones.
