.. _sec-software-overview:

****************************
ShakeMap Software Overview 
****************************

ShakeMap is a collection of modules written in PERL and C.  PERL is a powerful, freely 
available scripting language that runs on all computer platforms.  The collection of PERL 
modules allows the processing to flow in discrete steps that can be run collectively or 
individually.  Within the PERL scripts, other software packages are called, specifically 
packages that enable the graphics and much of the heavy grid-based computation.  For 
instance, maps are made using the Generic Mapping Tool (GMT; :ref:`Wessel and Smith, 
1991) <gmt_ref>`, and the Postscript output from GMT is converted to JPEG format using 
`ImageMagick <http://www.imagemagick.org/script/index.php>`_ or 
`GraphicsMagick <http://www.graphicsmagick.org/>`_.  In the design of ShakeMap, 
all components are built 
from freely available, open-source packages. 

While the PERL language is not the fastest possible way to implement ShakeMap, we 
note that much of the heavy computational load is handled by highly optimized programs 
(usually written in C) that are called from within the PERL programs. Even for networks 
with hundreds of stations over large regions, ShakeMap takes only a minute or so to run 
on a modern computer (and much of that time is spent in product generation, e.g., 
converting PostScript images to JPEG---something that would be very difficult to 
optimize further).

To enable customization for specific earthquakes or for different regions, each ShakeMap 
module has an accompanying collection of configuration files.  For example, in these 
files, one assigns the regional geographic boundaries and mapping characteristics to be 
used by GMT, which ground motion prediction equation (GMPE) to use, where and how 
to transfer the maps, email recipient lists, and so on.  Specific details about the software 
and configuration files are described in detail in the :ref:`Software Guide <sm35_software_guide>`.

With standardization in GIS and web application interfaces (API), several aspects of the 
ShakeMap system could be accomplished within GIS applications, but the open-source, 
freely available nature of GMT combined with PERL scripting tools allows for a flexible 
and readily available ShakeMap software package.  Nonetheless, we do generate a 
number of GIS product formats for that widespread user group as described in the :ref:`users-guide`.

