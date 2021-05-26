.. _software-guide-4:

###############
Software Guide 
###############

ShakeMap is a collection of modules written in Python and some C.  Python is a powerful, freely 
available scripting language that runs on all computer platforms.  The collection of Python 
modules allows the processing to flow in discrete steps that can be run collectively or 
individually.  Within the Python modules, other modules are imported to efficiently
perform many of the operations and to ease the development. An installation script
allows for the simple installation of Python and the ShakeMap dependencies in a
Conda virtual environment.  In the design of ShakeMap, all components are built 
from freely available, open-source packages. 

To enable customization for specific earthquakes or for different regions, ShakeMap 
collection of configuration files.  For example, in these 
files, one assigns the regional geographic boundaries and mapping characteristics,
which ground motion prediction equations (GMPEs) to use, where and how 
to transfer the maps, email recipient lists, and so on.  Specific details about the software 
and configuration files are described in detail in the
:ref:`Installation and Configuration <sec-installation-4>` section.

ShakeMap generates a wide range of products designed to suit most users'
needs.  See the section :ref:`sec-products-4` in the Users Guide
for a description of most of these products..

.. toctree::
   :maxdepth: 1

   sg_architecture.rst
   sg_installation.rst
   sg_queue.rst
   sg_input_formats.rst
   sg_output_formats.rst
   sg_operational.rst
   sg_contributing.rst
