.. _sec_disclaimers:

================
Disclaimers
================

General Disclaimer
---------------------------
.. warning:: Some USGS information accessed through this page may be preliminary in nature and presented prior to final review and approval by the Director of the USGS. This information is provided with the understanding that it is not guaranteed to be correct or complete, and conclusions drawn from such information are the sole responsibility of the user. In addition, ShakeMaps are automatic, computer-generated maps and have not necessarily been checked by human oversight, so they may contain errors. Further, the input data is raw and unchecked, and may contain errors.

* Contours can be misleading, since data gaps may exist. Caution should be used
  in deciding which features in the contour patterns are required by the data.
  Ground motions and intensities can vary greatly over small distances, so these
  maps are only approximate.

* Locations within the same intensity area will not necessarily experience the
  same level of damage, since damage depends heavily on the type of structure,
  the nature of the construction, and the details of the ground motion at that
  site. For these reasons, more or less damage than described in the intensity
  scale may occur. The ground motion levels and descriptions associated with
  each intensity value are based on recent damaging earthquakes. There may be
  revisions in these parameters as more data become available or due to further
  improvements in methodologies.

* Large earthquakes can generate very long-period ground motions that can cause
  damage at great distances from the epicenter; although the intensity estimated
  from the ground motions may be small, significant effects to large structures
  (e.g., bridges, tall buildings, storage tanks) may be notable.

* The utilization of DYFI data on ShakeMap, in addition to using recorded peak
  ground motions, is standard-operating-procedure for some ShakeMap
  operations, including the Global ShakeMap (GSM) system operated out of the
  USGS NEIC. The algorithmic strategy for including these data in ShakeMap is
  documented in the Technical Guide. As described by Wald et al.
  (2011), the ultimate discretion to use, filter, or overrule specific
  internet-based or historic intensities (or specific suspect strong motion
  data, for that matter) remains with the ShakeMap operators. A number of
  filtering and quality control strategies are in place (e.g., Wald et al.,
  2011), but erroneous or suspect data can not always be culled immediately.
  While we make efforts to provide consistent quality control of the data, the
  DYFI system depends upon open, citizen-science based input from the public.
  A number of studies have shown these data to be generally reliable, but the
  data reliability may vary from event to event. Moreover, macroseismic
  intensities are fundamentally non-unique: Differing polygonal aggregations
  for computing Community Decimal Intensities (CDI, using geocoded boxes or
  ZIP codes, for example, Wald et al., 2011) may yield varying intensity
  values at specific locations. [Historic or modern MMI or EMS-98 intensity
  assignments are also non-unique; the assignment can vary from expert to
  expert, the selection of areas may vary, and occasionally different
  structure types may indicated alternative intensities]. Changes to the
  size of the areas used to aggregate CDIs often trade off a greater number 
  of responses per polygon (hence greater confidence in the derived
  intensity) against a more precise spatial location. DYFI data are
  routinely used on the GSM systems and other regional ShakeMap systems of the
  ANSS. DYFI data are not currently (as of 2016) used in the Northern or 
  Southern California ShakeMap systems, in part due to the adequacy of
  strong-motion station coverage there.


ShakeMap Update Policy
---------------------------------------------------

.. warning:: ShakeMaps are preliminary in nature and will be updated as data arrives from a variety of distributed sources. Our practice is to improve the maps as soon as possible, but there are factors beyond our control that can result in delayed updates.

* For regions around the world, where there are insufficient near--real-time
  strong-motion seismic stations to generate an adequate strong-ground-motion
  data-controlled ShakeMap, we can still provide a very useful estimate of the
  shaking distribution using the ShakeMap software. Site amplification is
  approximated from a relationship developed between topographic gradient and
  shear-wave velocity. Additional constraints for these predictive maps come
  primarily from  additional earthquake source information, particularly fault
  rupture dimensions, observed macroseismic intensities (including via the USGS
  "Did You Feel It?" system), and observed strong ground motions, when and where
  available.
    
* There is no formal “final” version of any ShakeMap. Version Control is up to
  users. ShakeMap version numbers and timestamps are provided on the maps, 
  webpages, grid files, and metadata.

* Our strategy is to update ShakeMaps as warranted from a scientific
  perspective. We reserve the option to update ShakeMaps as needed to add data
  or to improve scientific merit and/or presentation of the maps in any way
  beneficial. This most typical update is after new data arrive, finite-fault
  models get established or revised, magnitude gets revised, or as improved site
  amplification maps, ground motion prediction equations, or even interpolation
  or other procedures become available. 

.. sidebar:: Recent ShakeMap update examples:

  * For the very deadly 2008 Wenchuan, China, earthquake, the Chinese strong-motion data were not made available for several months. 
  * For the 2011 Tohoku, Japan earthquake, the magnitude was updated from 7.9 to 8.9 over the course of the first hour after origin time. The Japanese strong-motion data processing center was impacted by the earthquake, yet they provided data for nearly a thousand seismic stations within several days after the earthquake. These vital data were added to the ShakeMap as soon as they became available.
  * Due to telemetry limitations, some important seismic station data for the 2014 American Canyon, California, earthquake came in minutes, hours, and as late as four days after the event. The data were added to the ShakeMap soon after they were received and processed. The magnitude also changed from an initial M5.7 to M6.0, and this, too, affected the ShakeMap. Lastly, the causative fault location was added by the Northern California ShakeMap operators several days after the earthquake.

**Updates to Online Maps**

   * **Real-time ShakeMap Updates**. Changes can be tracked with the ShakeMap
     version numbers and timestamps found in the metadata, the *info.xml* and
     *grid.xml* files, and on the maps themselves (time generated). The *info.xml*
     file contains timestamps, number of stations used, GMPE information, and
     many other attributes that could have changed from version to version.
     Often a text statement is provided that notes significant changes for a
     particular version. 

   * **ShakeMap Atlas Updates**. The ShakeMap Atlas uses version numbers
     for each Atlas event; the overall Atlas collections is also Versioned.
     Currently ShakeMap Atlas Version 2.0 is online in the ComCat database, and
     the older Atlas (Version 1.0) can be found online on the `legacy ShakeMap
     Archives pages <http://earthquake.usgs.gov/earthquakes/shakemap/>`_.

   * **Scenario Revisions**. ShakeMap Scenario collections uses version
     numbers for each event; the overall scenario collections is named
     according to their source. Currently
     ShakeMap Atlas Version 2.0 is online in the USGS `Comprehensive Catalogue
     (ComCat) Earthquake database
     <http://earthquake.usgs.gov/earthquakes/search/>`_. Some older
     scenarios are online on the `legacy ShakeMap Archives pages
     <http://earthquake.usgs.gov/earthquakes/shakemap/>`_. Scenario ShakeMaps
     will be revised when the underlying probabilistic seismic map fault
     segmentation or other particulars (like GMPE selection) change. Older
     versions will be archived in `ComCat
     <http://earthquake.usgs.gov/earthquakes/search/>`_.


	
