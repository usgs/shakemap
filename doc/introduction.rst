.. _introduction:

####################
Introduction
####################

.. only:: html

   This online ShakeMap Manual (:ref:`Worden and Wald, 2016 <worden2016b>`), is a 
   significant update of the
   original (:ref:`Wald et al., 2005 <wald2005>`) ShakeMap Manual. The
   Manual is now dynamic and served online. We employ Python document 
   generator `Sphinx <http://sphinx-doc.org/>`_, with the source
   under `GitHub <http://github.com/>`_ version control. These web pages
   are rendered with the `Sphinx RTD theme <http://github.com/snide/sphinx_rtd_theme>`_.
   A `PDF version <_static/ShakeMapManual.pdf>`_ of this manual is available.

.. only:: latex

   This ShakeMap Manual (:ref:`Worden and Wald, 2016 <worden2016b>`), is a 
   significant update of the
   original (:ref:`Wald et al., 2005 <wald2005>`) ShakeMap Manual.  We employ Python document 
   generator `Sphinx <http://sphinx-doc.org/>`_, with the source
   under `GitHub <http://github.com/>`_ version control.
   The online version of the manual is available at http://usgs.github.io/shakemap.

.. only:: html

   .. figure::  _static/Napa_ShakeMap_cover.*
      :scale: 50%
      :alt: ShakeMap Napa
      :align: right
      :target: Napa ShakeMap Example (URI or reference name)

      2014 M6.0 Napa, CA, earthquake intensity ShakeMap.

.. only:: latex

   .. figure::  _static/Napa_ShakeMap_cover.*
      :alt: ShakeMap Napa
      :scale: 70%
      :align: center
      :target: Napa ShakeMap Example (URI or reference name)

      2014 M6.0 Napa, CA, earthquake intensity ShakeMap.

`ShakeMap® <http://earthquake.usgs.gov/shakemap/>`_, 
developed by the U.S. Geological Survey (USGS), facilitates communication of 
earthquake information beyond just magnitude and location. By rapidly mapping out 
earthquake ground motions, ShakeMap portrays the distribution and severity of shaking. 
This information is critical for gauging the extent of the areas affected, determining which areas 
are potentially hardest hit, and allowing for rapid estimation of losses. Key to 
ShakeMap's success, algorithms were developed that take advantage of any high-quality 
recorded ground motions---and any available macroseismic intensity data---to provide 
ground-truth constraints on shaking. Yet ShakeMap also utilizes best practices
for both interpolating recordings and---critically---providing
event-specific estimates of shaking in areas where observations are sparse or nonexistent. Thus, ShakeMap portrays 
the best possible description of shaking by employing a combination of recorded and 
estimated shaking values. 

This Manual provides background on technical aspects of ShakeMap including: 1) information on 
the wide range of products and formats ShakeMap produces, 2) the uses of these products, 
and 3) guidance for 
ShakeMap developers and operators. Readers interested in understanding how 
ShakeMaps works can navigate to the :ref:`technical-guide`. Those who want to use 
ShakeMap products and understand their varied forms can jump to the :ref:`users-guide`. 
If your goal is to install and operate ShakeMap, see the :ref:`software-guide`. The
:ref:`software-guide` also points users to the ShakeMap software distribution and 
`Software Guide <http://usgs.github.io/shakemap/_static/SoftwareGuideV3_5.pdf>`_.
