.. _sec-tg-parameters:

*************************************
Recorded Ground-motion Parameters
*************************************

Data Acquisition
=================

ShakeMap requires estimates of magnitude, location, and (optionally, but preferably) 
shaking IMs at seismic stations. As such, ShakeMap has been interfaced with 
several types of seismic processing systems in wide use at numerous networks across the U.S. and 
around the world, including `Antelope <http://antelopeusersgroup.org/>`_, 
`SeisComP3 <https://www.seiscomp3.org/>`_, and `AQMS <http://www.isti.com/products/aqms/>`_. 
The ShakeMap system, 
however, is a stand-alone software package and is really a passive consumer of seismic 
data. In other words, the ShakeMap software itself contains no data acquisition component. It is 
assumed that ShakeMap earthquake source information and parametric data are packaged 
for delivery to ShakeMap and that that delivery will trigger a ShakeMap run. The 
required format is an XML format, as fully described in the :ref:`Software Guide <sm35_software_guide>`. 
Some programs are provided to convert ASCII text and other formats to the required input
XML.
It is assumed that station data delivered to ShakeMap are free-field sites that have been 
vetted by the contributing network. Each station must have stand-alone metadata 
describing its station location, contributing network, channel, and location code. While some 
additional outlier and data quality checks are performed within ShakeMap (see 
:ref:`sec_shakemap_processing`), it is assumed that this is primarily the 
responsibility of the contributing seismic network.

For global and historic earthquake ShakeMap generation, we have developed scripts to 
preprocess various forms of seismic waveform (as well as macroseismic) data which are 
openly available around the world. For example, we provide a Python 
script `getstrong.py <https://github.com/mhearne-usgs/smtools>`_
that runs independently of ShakeMap, as described in the :ref:`software-guide`.

For illustrative purposes, we describe the data acquisition for the seismic system in 
Southern California, a component of the California Integrated Seismic Network (`CISN <http://www.cisn.org>`_). 
For perspective, as of 2015, there were nearly 800 real-time stations jointly operated with 
a collaboration between the USGS and the California Institute of Technology (Caltech). In addition, the California Geological Survey (CGS) 
contributes nearly 350 strong-motion stations in near real-time, utilizing an automated 
telephone dial-up procedure (:ref:`Shakal et al, 1998 <shakal1998>`), and the USGS National Strong Motion 
Instrumentation Program (NSMP) contributes dial-up station parameters as well, with 
nearly 50 stations in Southern California alone.  Lastly, the
`"NetQuakes" <http://earthquake.usgs.gov/monitoring/netquakes/map>`_ program, a 
relatively low-cost seismograph that the USGS installs in homes, businesses, buildings, and 
schools, contributes close to 100 additional stations in Southern California. 

Generation of ShakeMap in Southern California is automatic, triggered by the event 
associator of the seismic network.  Within the first two minutes of an earthquake, ground-motion 
parameters are available from the USGS-Caltech component of the network, and 
within several minutes most of the important near-source CGS and NSMP stations contribute; a 
more complete CGS and NSMP contribution is available within approximately 
ten minutes of an event.  Initial maps are made with the real-time component of CISN as 
well as any of the available dial-sites, and they are updated automatically as more data 
are acquired.

Derived Parametric Ground-motion Values
=========================================

Parametric data from stations serving ShakeMap should include peak ground acceleration 
(PGA), peak ground velocity (PGV), and peak response spectral acceleration amplitudes 
(at 0.3, 1.0, and 3.0 sec).  Often, parametric values are derived continuously, using 
recursive time-domain filtering as described by :ref:`Kanamori et al. \(1999\) <kanamori1999>`.  Otherwise 
parameters are derived from post-processing as described by :ref:`Shakal et al. \(1998\) <shakal1998>` and 
:ref:`Converse and Brady \(1992\) <converse1992>`.

ShakeMap will run successfully with no (or limited) parametric data, for example if only 
PGA values are available at each station. Default logic is employed to provide reasonable 
behavior for estimating intensities from PGA alone, bias correction, and interpolation 
(see following sections). Likewise, for smaller-magnitude earthquakes, spectral values 
can be noisy, so operators often omit the generation spectral maps below a lower 
magnitude threshold (about M4); this can be done with simple command-line options. 

For all maps and products, the motions depicted are peak values as observed; that is, the 
maximum value observed on the two horizontal components of motion.  Many engineers 
are accustomed to analyses employing the geometric mean of the horizontal peak-ground 
motions, but that parameter is not computed by ShakeMap. More description and justification 
for this strategy is given in the section :ref:`peak-values-vs-mean`. It should be noted, 
however, that conversions from peak to geometric mean [or orientation-independent 
parameterizations (:ref:`Boore, 2010 <boore2010>`)] are available 
(e.g., :ref:`Beyer and Bommer, 2006 <beyer2006>`).

Macroseismic Intensity
===========================

ShakeMap also (optionally) accepts input data in the form of observed macroseismic 
intensity (MMI, MCS, etc.). As with peak ground motion parameters from seismic 
stations, ShakeMap expects specific file formats (XML) and site metadata for 
macroseismic data (see the :ref:`Software Guide <sm35_software_guide>`).

Intensity data can fill important gaps where ground-motion recordings are not available, 
and often provide the only control in sparsely instrumented areas. This is particularly true 
for historic earthquakes, for which macroseismic data provide important constraints on 
shaking intensities. As later discussed, the ShakeMap Atlas (:ref:`Allen et al., 2008 <allen2008>`, :ref:`2009a <allen2009a>`; 
:ref:`Garcia et al., 2012a <garcia2012a>`) is a collection of important historic earthquake shaking maps which 
are now widely used for scientific analyses and for loss model calibration (e.g., :ref:`Wald et 
al., 2008 <wald2008>`; :ref:`Jaiswal and Wald, 2010 <jaiswal2010>`; :ref:`Pomonis and So, 2011) <pomonis2011>`. 

The most common source for immediate post-earthquake intensity data is the USGS's 
“Did You Feel It?” (DYFI) system (:ref:`Wald et al., 2011 <wald2011c>`), though similar systems are 
available in several countries. However, traditionally assigned intensities may be used as 
well. DYFI data can be programmatically retrieved from the USGS's database and 
formatted for ShakeMap input using the ShakeMap program *getdyfi*, making it 
especially easy to incorporate into the ShakeMap data input stream. 

Macroseismic intensity data can also be an important constraint on peak ground motions, 
since ground motion amplitudes can be derived from intensity through the use of a suitable Ground-Motion/Intensity 
Conversion Equation (GMICE). Because a GMICE represents a statistical (probabilistic) 
relationship, the conversion to and from intensity has a higher uncertainty than direct 
ground-motion observation. ShakeMap accounts for this higher uncertainty by down-
weighting converted observations in the interpolation process, as discussed in the 
:ref:`sec_interpolation` section.

A variety of GMICEs are available with the ShakeMap software distribution, both for 
MMI---based on :ref:`Wald, et al. \(1999b\) <wald1999b>`, :ref:`Worden, et al. \(2012\) <worden2012>`, 
and :ref:`Atkinson and Kaka \(2007\) <atkinson2007>`, among others---and for MCS---based on :ref:`Faenza and Michilini \(2010\) <faenza2010>`. Operators are 
encouraged to explore the need to develop their own relationships based on data covering 
their own operational area as GMICEs have been shown to have regional dependencies 
(e.g., :ref:`Caprio et al., 2015 <caprio2015>`). A complete list of GMICEs currently employed by ShakeMap is 
provided in the :ref:`Software Guide <sm35_software_guide>`.

We have implemented a convention for maps and regression plots that seismic stations 
are represented with triangles and macroseismic data are depicted with circles (see :num:`Figure 
#figure1-1`, for example). This convention is forward-looking: not all seismic networks were 
currently following this convention at the time of this writing.  

.. _figure1-1:
 
.. figure:: _static/Figure_1_1.*
   :align: left
   :width: 650px

   Intensity ShakeMap from the 2014 M6.0 American Canyon (Napa Valley), 
   CA earthquake. Strong motion data (triangles) and intensity data (circles) are color-coded 
   according to their intensity value, either as observed (for macroseismic data) or as converted 
   by :ref:`Wald et al. \(1999b\) <wald1999b>` as shown in the legend. The north-south black line indicates the 
   fault location, which nucleated near the epicenter (red star). Note: Map Version Number reflects 
   separate offline processing for this Manual. 

:num:`Figure #figure-hawaii-interactive` shows a different representation of the 
intensity map  on the newer, interactive maps on the USGS web site.

.. _figure-hawaii-interactive:

.. figure:: _static/Hawaii_interactive.*
   :align: left
   :width: 650px
   :scale: 90 %
  
   Intensity ShakeMap from the 2006 M6.7 Kahola Bay, HI earthquake. 
   Contours indicate intensities; strong motion data (triangles) and intensity data (circles) are color-
   coded according to their intensity value, either as observed (for macroseismic data) or as 
   converted by :ref:`Worden et al. \(2012\) <worden2012>`. Inset on lower map shows 
   pop-up station information. 

