.. _sec-philosophy:

**********************************************************
Philosophy of Estimating and Interpolating Ground Motions
**********************************************************

The overall strategy for the deployment of stations under the ANSS implementation plan 
relies on dense instrumentation concentrated in urban areas with high seismic hazards 
(:ref:`USGS, 1999 <usgs1999>`) and fewer stations in outlying areas.  Based on this philosophy, and when 
fully deployed, maps generated in these urban regions are expected to be most accurate 
where the population at risk is the greatest, and therefore, where emergency response and 
recovery efforts will likely be most urgent and complex.  

Even so, significant gaps in the observed shaking distribution will likely remain, 
especially in the transition from urban to more rural environments. Likewise, many 
critical facilities and lifelines are widely distributed, away from population centers and 
their dense seismic sensor networks.  Thus, as a fundamental strategy for ShakeMap, we 
have developed algorithms to best describe the shaking in more remote areas by utilizing 
a variety of seismological tools.  In addition to the areas without sufficient 
instrumentation where we would like to estimate motions to help assess the situation, and 
as a fail-safe backup, it is also useful to have these algorithms in place in 
the event of potential communication dropout from a portion of the network.  The same 
tools are, in fact, beneficial for interpolating between observations (i.e., seismic stations) even 
within densely instrumented portions of the networks.

If there were stations at each of the tens of thousands of map grid points needed to 
adequately portray shaking, then the creation of shaking maps would be relatively simple.  
Of course, stations are not available for the overwhelming majority of these grid points, and in many cases grid 
points may be tens of kilometers or more from the nearest reporting station.  The overall mapping 
philosophy is then to combine information from individual stations, site amplification 
characteristics, and ground-motion prediction equations for the distance to the hypocenter 
(or to the causative fault) to create the best composite map.  The procedure should 
produce reasonable estimates at grid points located far from available data while 
preserving the detailed shaking information available for regions where there are stations 
nearby.

It should be mentioned that mathematically, or algorithmically, geospatial interpolation 
can take many forms. There are some good reasons to employ geospatial kriging-with-a-trend. 
However, the complexity of the trends (GMPE, as well as inter-event bias 
corrections per Intensity Measure or IM), the use of multiply-weighted strong-motion and macroseimic 
data, and the real-time nature of the processing require other considerations. Effectively, 
the approach ShakeMap currently employs for interpolation (:ref:`Worden et al., 2010 <worden2010>`), which 
employs a predetermined spatial correlation function, is broadly analogous to `kriging-with-a-trend <https://en.wikipedia.org/wiki/Kriging>`_
mathematically. We address this possibility further in :ref:`future-directions`.

Estimating motions where there are few stations, and then interpolating the recordings and 
estimates to a fine grid for mapping and contouring, requires several steps. In the 
following sections, we describe the process from input to final interpolated grid. Where 
beneficial, we illustrate the effects of key steps with example ShakeMap figures.

