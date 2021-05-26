.. _sec_disclaimers-4:

================
Disclaimers
================

General Disclaimer
---------------------------
.. warning:: Some USGS information accessed through this page may be
             preliminary in nature and presented prior to final review and
             approval by the Director of the USGS. This information is
             provided with the understanding that it is not guaranteed to
             be correct or complete, and conclusions drawn from such
             information are the sole responsibility of the user. In
             addition, ShakeMaps are automatic, computer-generated maps and
             have not necessarily been checked by human oversight, so they
             may contain errors. Further, the input data is raw and
             unchecked, and may contain errors.

* Contours can be misleading, since data gaps may exist. Caution should be
  used
  in deciding which features in the contour patterns are required by the data.
  Ground motions and intensities can vary greatly over small distances, so
  these
  maps are only approximate. In addition, the contours provided in the JSON
  files have been smoothed, and are thus not suitable for anything but 
  plotting and approximate analysis.

* Locations within the same intensity area will not necessarily experience the
  same level of damage, since damage depends heavily on the type of structure,
  the nature of the construction, and the details of the ground motion at that
  site. For these reasons, more or less damage than described in the intensity
  scale may occur. The ground motion levels and descriptions associated with
  each intensity value are based on recent damaging earthquakes. There may be
  revisions in these parameters as more data become available or due to
  further improvements in methodologies.

* Large earthquakes can generate very long-period ground motions that can
  cause damage at great distances from the epicenter; although the intensity
  estimated from the ground motions may be small, significant effects to
  large structures (e.g., bridges, tall buildings, storage tanks) may be
  notable.

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

.. warning:: ShakeMaps are preliminary in nature and will be updated as
             data arrives from a variety of distributed sources. Our
             practice is to improve the maps as soon as possible, but there
             are factors beyond our control that can result in delayed
             updates.

Also see: :ref:`ug-sec-shakemap-operations-4`, in particular
:ref:`subsec-sources-versions`.
	
