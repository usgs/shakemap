.. _future-directions:

####################
Future Directions
####################

ShakeMap is a continual work-in-progress. We note several ongoing developments and "To-Do" lists. Please make suggestions if you would like to weigh in.


.. sidebar:: Feedback? Suggestions?
	     
    Please `let us know <http://earthquake.usgs.gov/contactus/?to=cbworden,klin,wald&subject=ShakeMap>`_ what you're thinking!

Research & Development
---------------------------------
**Feature Requests:**

* Add NGA-West2, NGA-East, and NGA-Subduction GMPEs, including basin terms for NGA-West 2 GMPEs.
* Improved and additional site amplification approaches and tables, in addition to native GMPE (Vs30) site corrections (e.g., Choi and Stewart,
  2005).
* R&D to improve PGV-to-MMI conversion for large-magnitude and high-velocity ranges. May require switch to converting long-period spectral acceleration
  to MMI. Simulated ground motion time histories may be useful to augment sparse data at high PGV/MMI. 
* Consideration of vector-component IMs, static displacements, and duration-based IMs (Arias Intensity; Cumulative Average Velocity, or CAV) [#f1]_

**In Progress:**

* Spatial variability. Implement optimization methods to compute the spatial correlation field for ShakeMap
  using successive conditional simulations (`Verros et al., 2016) <verros2016>`_.
  Deliver ShakeMap scenarios with multiple realizations of variability.  
* Directivity. Update `Rowshandle (2010) <rowshandle2010>`_ model and implement selected NGA-West2 models.
* Landslide and liquefaction likelihood grid (*sechaz.grid.xml*). Computing probability and distribution of landsliding and liquefaction per
  ShakeMap grid cell. Delivery via Product Distribution Layer (PDL) for ShakeCast, PAGER, and open access.
* Scenario update for all U.S. regions. Delivery to ComCat/Web and allow users a variety of search capabilities (site- or fault-specific).
* Interactive (dynamic) webpages plots (regression, bias, outliers, station amplitudes).
* Improved content and rendering of ShakeMap metadata (*info.xml*; see :ref:`Thompson et al., 2016 <thompson2016>`).  


Software: ShakeMap 4.0 (Python)
-----------------------------------
The release of ShakeMap version 4.0 will represent a significant departure from previous versions. All of the important computational modules have been refactored into the Python programming language and make use of the tools in the widely available Python “scientific distribution”. The core ShakeMap code, approaching fifteen years old, was due for a major overhaul---to more organically incorporate the many extensions that had been added over its lifetime, and to facilitate several new demands from ShakeMap software and ShakeMap’s expanded role as a global provider of post-earthquake information, earthquake scenarios, and inputs for loss-modeling software.  

One of the advantages of the rewrite of ShakeMap into the Python language was the availability of the GEM OpenQuake (OQ) library of Ground Motion Prediction Equations (GMPEs). The OQ hazard library provided us with a broad range of well-tested, high-performance, open-source GMPEs. Due to constraints imposed by the software architecture of earlier implementations of ShakeMap, the development of GMPE modules was time-consuming and difficult, which restricted the number and timeliness of the available modules. The OQ library provides a broad array of current GMPE modules, as well as a framework for easily adding new modules (whether by GEM or ShakeMap staff), jumpstarting our efforts to re-implement ShakeMap.

The OQ hazard library also provides supporting functions for using the GMPE modules, including a set of software classes for computing the various distance measures required by the GMPEs. The ShakeMap fault model, however, was somewhat more general than allowed for by the OQ planar surface modules, so we sub-classed the OQ “surface” class and implemented our own high-performance module. The open-source, cooperative nature of the OQ project allowed us to contribute our new module back to the OQ repository, and thus make it available to other users.

In addition to making use of the GEM OQ library, there are a number of other advantages to using Python for an application like ShakeMap.  The dynamic nature of the language means that development time is much reduced, allowing a small team to generate useful code in a short amount of time.  Also, there is an active scientific computing Python community that has created many tools that solve common problems, including an array object for vectorized operations, input/output routines for common data formats, and plotting/mapping libraries.  These tools further help to reduce development time and effort.
[**Delivery Date: 2016**] 

.. [#f1] We are continuously considering the use of additional ground-motion parameters (IMs)
	 for ShakeMap. However, any such additions cannot be made lightly. In part, this is
	 due to the fact that this requires upgrading process seismic network processing streams
	 that produce parametric and these processes vary significantly among ANSS data sources.    

