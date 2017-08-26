.. _sec-architecute-4:

****************************
Architecture
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

.. figure:: _static/consumer.*
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

