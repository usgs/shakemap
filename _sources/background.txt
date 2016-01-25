.. _sec_background:

===================
Background
===================

Until the development of ShakeMap, the most common information available
immediately following a significant earthquake was typically its magnitude and
epicenter.  However, the damage pattern is not a simple function of these two
parameters alone, and more-detailed information must be provided to properly
assess the situation.  For example, for the 1971 M6.7 San Fernando,
California earthquake, the northern San Fernando Valley was the region with the
most damage, even though it was more than 15km from the epicenter.  Likewise,
areas strongly affected by the 1989 M6.9 Loma Prieta and 1994 M6.7 Northridge, California,
earthquakes that were either distant from
the epicentral region or out of the immediate media limelight were not fully
appreciated until long after the initial reports of damage. The full extent of
damage from the 1995 M6.9 Kobe, Japan, earthquake was not recognized by
the central government in Tokyo until many hours later (e.g., :ref:`Yamakawa, 1998 <yamakawa1998>`),
seriously delaying rescue and recovery efforts.

In contrast, a ShakeMap is a representation of actual
ground shaking produced by an earthquake. The
information it presents is different from the earthquake magnitude and epicenter
that are released after an earthquake, because ShakeMap focuses on the
ground shaking produced by the earthquake, rather than the parameters describing
the earthquake starting point (its hypocenter) and size (magnitude). So,
although an earthquake has one magnitude and one epicenter, it produces a range
of ground shaking levels at sites throughout the region, depending on distance
from the earthquake fault that ruptured, the rock and soil conditions at sites,
and variations in the propagation of seismic waves from the earthquake due to
complexities in the structure of the Earth's crust. 

Part of the strategy for generating rapid-response ground motion maps was to
determine the best format for reliable presentation of the maps given the
diverse audience, which includes scientists, businesses, emergency response
agencies, media, and the general public.  In an effort to simplify and maximize
the flow of information to general users, we have developed a means of generating
not only PGA and PGV maps, but also an
instrumentally derived estimated Modified Mercalli Intensity (MMI) map.  This
“instrumental intensity” map makes it easier to relate the recorded
ground motions to the expected felt area and damaging shaking distribution. At the same time,
we preserve a full range of utilities of recorded ground-motion data by
producing maps of response spectral acceleration, which are not particularly
useful to the general public, yet which provide fundamental data for loss
estimation and earthquake engineering assessments.

As mentioned, ShakeMap provides maps of **peak** ground-acceleration, velocity, and spectral
acceleration as well as MMI. Intensity ShakeMaps
depict estimated intensities derived from peak ground motions as well as
(optionally) from reported intensities. Intensity maps make it easier to relate
the recorded and estimated ground motions to the expected felt and damage
distributions. Intensities are estimated from ground shaking 
using equations based on analyses
of intensities reported near recorded seismic stations for past
earthquakes, e.g., :ref:`Wald et al. \(1999b\) <wald1999b>` or
:ref:`Worden et al. \(2012\) <worden2012>`. The
legend on the ShakeMap indicates which relationship was used to estimate
intensities from ground motions and vice versa (see the ShakeMap
:ref:`technical-guide` for more details).

Station locations are the best indicator of where the map is most accurate: near
seismic stations, the shaking is well constrained by data; far from such
stations, the shaking is estimated using standard seismological inferences and
geospatial interpolation. Details about the interpolation; uncertainty maps; and
codes for the seismic station components, network contributors, and potential
outlier or clipped flag codes are provided in the :ref:`technical-guide`. Peak
horizontal acceleration and spectral acceleration values are in units of
percent-g (or %g, where g = acceleration due to the force of gravity = 981cm/s/s). The
peak values of the vertical components are not used in the construction of the
maps. Peak velocity values are given (in cm/s) at each station. Acceleration
spectra are the response of a 5% critically damped, single-degree-of-freedom
oscillator to the recorded ground motion at three reference periods: 0.3, 1.0,
and 3.0 sec. 
