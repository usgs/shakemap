.. _tg-sec-choice-of-parameters:

****************************************************************
Discussion of Chosen Map Parameters (Intensity Measures)
****************************************************************

.. _peak-values-vs-mean:

Use of Peak Values Rather than Mean
============================================

With ShakeMap, we chose to represent peak ground motions as recorded. We depict the 
larger of the two horizontal components, rather than as either a vector sum, or as a 
geometric mean value. The initial choice of peak values was necessitated by the fact that 
roughly two-thirds of the TriNet (now the Southern California portion of CISN) strong-
motion data (the California Geological Survey, or CGS, data)
are delivered as peak values for individual components of 
motion, that is, as parametric data rather than waveforms. This left two options: providing peak 
values or median of the peak values---determining vector sums of the two horizontal 
components was not an option, because the peak values on each component do not 
necessarily occur at the same time.  A useful strategy going forward may be to employ 
the 50th percentile of the response spectra over all non-redundant rotation angles 
(RotD50; :ref:`Boore, 2010 <boore2010>`), which is becoming a standard
for "next-generation" 
GMPEs (:ref:`Abrahamson et al., 2014 <abrahamson2014>`), or on another
agreed-upon vector-component 
calculation. (See :ref:`future-directions`). However, such changes would
require adoption of these 
calculations on the part of each contributing seismic network, and would thus require 
consensus (and substantial software upgrades) all around. 

Despite the common use of mean values in attenuation relations and loss estimation, we 
decided that computing and depicting median values, which effectively reduces 
information and discards the largest values of shaking, was not acceptable. This is 
particularly true for highly directional near-fault pulse-like ground motions, for which 
peak velocities can be large on one component and small on the other.  Mean values for 
such motions (particularly when determined in logarithmic space) can seriously underrepresent 
the largest motion that a building may have experienced, and these pulse-like motions are 
typically associated with the regions of greatest damage. Thus, we chose peak ground 
motions as the parameters to be mapped.

:ref:`Beyer and Bommer \(2006\) <beyer2006>` provide statistical relationships
to convert among median and 
peak parameters and between aleatory variability for different definitions of the 
horizontal component of motion. This is useful when approximating alternative 
components than those presented, but one must recognize that for any individual record, 
these statistics may or may not be representative. 

Initially, our use of PGA and PGV for estimating intensities was also simply practical. 
We were retrieving only peak values from a large subset of the network, so it was 
impractical to compute more specific ground-motion parameters, such as average-
response spectral values, kinetic energy, cumulative absolute velocities (CAV, :ref:`EPRI, 
1991 <epri1991>`), or the JMA intensity algorithm (:ref:`JMA, 1996
<jma1996>`). However, because
near-source strong ground motions are often dominated by short-duration, pulse-like
ground motions (usually associated with source directivity), PGV appears to be a robust measure 
of intensity for strong shaking. In other words, the kinetic energy (proportional to 
velocity squared) available for damage is well characterized by PGV. In addition, the 
close correspondence of the JMA intensities and peak ground velocity 
indicates that our use of peak ground velocities for higher intensities was 
consistent with the algorithm used by JMA. Work by :ref:`Wu et al. \(2003\)
<wu2003>` indicates a very 
good correspondence between PGV and damage for data collected on the island of Taiwan, 
which included high-quality loss data and densely sampled strong-motion observations 
for the 1999 Chi-Chi earthquake. More recent work on Ground-Motion/Intensity 
Conversion Equations (GMICEs) (e.g., :ref:`Atkinson and Kaka, 2007
<atkinson2007>`; :ref:`Worden, et al., 2012 <worden2012>`) has also 
confirmed the strong relationship between PGV and intensity. 

Nonetheless, for large, distant earthquakes, peak motions may be less informative, and 
spectral content and perhaps duration become key parameters.  Although we may eventually 
adopt corrections for these situations, it is difficult to assign intensities in such cases. For 
instance, it is difficult to assign the intensity in the zone of Mexico City where numerous 
high-rises collapsed during the 1985 Michoacan earthquake. There was obviously high-
intensity shaking for high-rise buildings; however, most smaller buildings were 
unaffected, suggesting a much lower intensity.  Whereas PGVs were 
moderate and would imply intensity VIII, resonance and duration conspired to cause a 
more substantial damage than peak values would indicate. Although this is, in part, a 
shortcoming of using peak parameters alone, it is more a limitation imposed by 
simplifying the complexity of ground motions into a single parameter. Therefore, in 
addition to providing peak ground-motion values and intensity, we are also producing 
spectral response maps (for 0.3, 1.0, and 3.0 sec). Users who can process this information 
for loss estimation will have a clearer picture than can be provided with maps of PGA 
and PGV alone. However, as discussed earlier, a simple intensity map is extremely useful 
for the overwhelming majority of users, which includes the general public and many 
people involved with the initial emergency response. 

