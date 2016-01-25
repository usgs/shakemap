.. _software-guide:

##########################################
Software & Implementation Guide 
##########################################

Installing and operating a ShakeMap system is a non-trivial endeavor. While the software 
can be easily obtained, and installed and configured with a few hours work, there are a 
great many issues that need to be addressed when an organization considers deploying a 
ShakeMap system. In this section, we will discuss a number of these issues.

***********************************
Implementation Considerations
***********************************

Seismic Network
=========================

The single most important prerequisite for operating a ShakeMap system is the existence 
of a seismic network capable of determining the location and magnitude of significant 
earthquakes within minutes of their occurrence. This network must also produce 
parametric data for ShakeMap within minutes of the earthquake. The parametric data 
should come from free-field strong motion sensors with real-time or near--real-time 
telemetry, and minimally consist of the horizontal components of Peak Ground 
Acceleration (PGA) and Peak Ground Velocity (PGV), as well as (ideally) 5% damped pseudo-
spectral acceleration (PSA) at 0.3, 1.0, and 3.0 sec periods.

In the absence of such a network, there is little reason to operate ShakeMap locally. The 
USGS operates a ShakeMap system that produces maps for significant global 
earthquakes within a few minutes of their occurrence, and incorporates any available data 
from the Global Seismic Network (GSN), “Did You Feel It?” (DYFI), and any other parametric data to which 
we have access. These maps are immediately delivered to the PAGER system and 
ShakeCast. It would be far simpler in many cases to obtain the 
ShakeMap products from the USGS (through the website or ShakeCast) and use them for 
local or regional response and recovery efforts, perhaps even hosting them on local 
servers, than to attempt to duplicate our existing system.

Even where a seismic network is operational, it may be preferable to work with the 
USGS to produce ShakeMaps. As the following sections will illustrate, running a 
ShakeMap system is a considerable operational undertaking, and in many instances it will 
strain the resources of local or regional networks. The USGS operates a fully-staffed, 
24/7 operations center for global earthquakes and produces ShakeMaps on a daily basis. 
If a regional network's parametric data can be made available to the USGS, the USGS's 
global ShakeMap system will include that data in the processing, and finished ShakeMap 
products can be delivered via ShakeCast for local consumption. Please contact us if you 
are interested in this kind of arrangement.

Seismic Stations and Parametric Data
===========================================

As mentioned above, seismic waveform or parametric data must be communicated to the 
processing center in real-time or near---real--time. The sensors must be "free field" (i.e., not 
located in, or adjacent to, large structures) and installed in accordance with modern standard 
practices. Strong-motion instruments are preferred for ShakeMap, as broadband 
instruments typically clip when ground motions reach the level of interest to ShakeMap. 
Co-located strong-motion and broadband instruments can, however, increase the dynamic 
range of a station and its overall usefulness, but the network operator is responsible for 
specifying the crossover from broadband to strong motion, and ShakeMap should be 
presented only with the amplitudes from the favored instrument. Off-scale, clipped, 
below-noise, or otherwise unreliable data must be flagged or omitted from the ShakeMap 
input files. Horizontal components are mandatory---ShakeMap does not use vertical 
components.

At a minimum, the network must produce PGA and PGV for the ShakeMap stations, but 
PSA (at 0.3, 1.0, and 3.0 sec) is also desirable. The algorithms 
used to compute the parametric data should be should be verified against known standards. 
ShakeMap can be quite helpful in 
highlighting grossly mis-calibrated stations, but it is best to find these (and more subtle) 
errors before a large earthquake strikes.

Intensity data from the USGS's DYFI system (or an equivalent 
regional internet intensity system) is acceptable as input to ShakeMap; however, care 
should be exercised with international DYFI data, as it is often aggregated by city and 
therefore may be too coarse-grained for large-scale maps.

Triggering ShakeMap
=========================

The operator must give careful consideration to the way ShakeMap will be triggered. 
Unless an operator is constantly standing by to run ShakeMap, the system must be 
automatic to be useful. The operator must select the boundaries of the region ShakeMap 
will cover and the minimum earthquake magnitude that will trigger a run. These choices 
will be influenced by the reliability of the network's earthquake location and magnitude 
data over the region in question. The operator must also consider that a large earthquake may 
occur outside their region of responsibility, but have effects inside their region.

The ShakeMap software is distributed with a program for responding to event triggers 
and queuing events to be processed, but it is only applicable to the AQMS network 
software. While this software can serve as a guide, users of other systems (e.g., Antelope, 
Earthworm) will need to develop their own triggering systems. 

Once triggered, ShakeMap will expect to find event and parametric data files in the 
event's input directory. Whether the amplitudes are retrieved from a database or stored in 
files, it is the operator's responsibility to provide ShakeMap with properly formatted 
files. Again, example programs are distributed with ShakeMap, but operators should 
anticipate that some coding may be necessary.

System Configuration
===========================

While it is relatively easy to install and run the ShakeMap software, a great deal of 
consideration is needed to properly configure the system for a particular region. 
ShakeMap presents the operator with a large number of configuration options, and these 
options will directly affect the accuracy and reliability of the products. The 
operator must select the proper GMPE or GMPEs for the region, and specify under which 
conditions they will be used. An appropriate Ground-Motion--Intensity Conversion 
Equation (GMICE) must likewise be chosen. The operator may also elect to use an Intensity 
Prediction Equation (IPE) or to use the default virtual IPE. For best results, the operator 
will need to provide a Vs30 grid, and while the `USGS's Vs30 server <http://earthquake.usgs.gov/hazards/apps/vs30/>`_ can supply such a 
grid based on topographic slope, a grid based on regional geology is preferable [#]_. The 
user must also decide whether to use GMPE-native or Borcherdt-style site correction 
factors. This is just a partial list of configuration choices. There are numerous other 
configurable parameters, all of which must be carefully chosen.

While the default ShakeMap webpages are configurable, they are fairly rudimentary and 
outdated. In addition, the webpages and ShakeMap products are all in English. Operators 
wishing more sophisticated webpages or non-English support can anticipate a substantial 
investment in bringing the system online. Some modifications of the ShakeMap software 
may be required for languages other than English, but such modifications may make it 
more difficult to update the software.

For further information on ShakeMap configuration, see the :ref:`Software Guide <sm35_software_guide>` and the 
ShakeMap configuration files themselves. 

Operations
=====================

Once an organization begins producing and distributing ShakeMaps, end users will begin 
to depend upon them and develop systems that incorporate ShakeMap products into their 
response and recovery operations. This means that it essentially becomes mandatory to 
produce ShakeMaps following significant earthquakes, and failure to do so is 
extremely conspicuous. A robust ShakeMap operation requires the 24/7 availability of 
operational personnel. If the facility is not continuously manned, on-call staff must be 
designated, and those staff must have the ability to access and operate the ShakeMap 
system remotely. Significant earthquakes almost always require some manual 
intervention (changing map scale, re-centering, addition of finite fault, 
inclusion/exclusion of data, etc.), and experienced personnel are required to evaluate the 
situation and perform the necessary tasks.

There are additional, more routine, operational considerations. An experienced 
seismologist should routinely review all of the ShakeMaps produced by the system and 
take action to correct any deficiencies. A network seismologist should also review the 
inputs and outputs of ShakeMap to insure that all stations are producing appropriate data. 
A ShakeMap operator should routinely review all ShakeMap processes, logs, databases, 
and outputs to insure the system is operating as expected. 

The ShakeMap software is usually updated a few times a year. These updates contain 
important bug fixes, new functionality, new products, and general improvements. An 
operator must review the change logs, decide when to apply the updates, and test the 
updated software before it is put into production mode. Occasionally it may be desirable 
to rerun earlier events or scenarios to take advantage of the capabilities of the new code. 

Hardware and software systems will need to be monitored and maintained for around-the-
clock availability. This includes not just the seismic network and ShakeMap systems, but 
also web servers and other network hardware and software required for delivering 
products to end users. The personnel responsible for these systems must be on-call and 
able to access the necessary systems remotely. Automatic monitoring of mission-critical 
hardware and software is strongly encouraged. These systems should also have several 
hours of backup power in case of an outage. Periodic outage tests should be 
conducted to ensure that all necessary systems remain operational.

As mentioned above, users can be expected to make use of ShakeMaps in a variety of 
ways. However, many organizations that could make use of ShakeMap products are 
unaware of ShakeMap and the ways it could serve their earthquake response and 
recovery needs. We have found that a sustained outreach effort is necessary to maximize 
the adoption of ShakeMap and, thus, its value to society. Potential end users include 
public utilities, government and private transportation companies, police and fire 
departments, regional and national emergency response organizations, private companies 
with distributed facilities (e.g., banks, chain stores, telecoms), insurance companies, 
investment houses, and many others. Not only can ShakeMap-improved response efforts 
benefit post-earthquake recovery, these organizations can provide much-needed support 
for network and ShakeMap operations. It is highly recommended that regional networks 
considering the implementation of ShakeMap develop a detailed outreach plan.

Scenarios
============

One important use of ShakeMap is the generation of earthquake scenarios. Scenarios are 
predictive maps of the potential shaking resulting from hypothetical future (or past) earthquakes. 
Scenarios can be used for planning exercises, public information, or research. Some 
users may request specific scenarios, but it is generally worthwhile to develop a suite of 
scenarios to cover the likely earthquake hazards of a region. At the USGS, we have begun 
using disaggregated hazard maps as the basis for our nationwide scenario project. In other 
words, we separate out the individual earthquakes (and causative faults) that together 
comprise the hazard in a probabilistic hazard map. The disaggregated maps represent the 
best scientific consensus of the probable earthquakes in a region, and should be sufficient 
for most uses. Requests for custom scenarios should be carefully evaluated. The 
earthquakes represented should be credible in terms of both the causative fault and the 
magnitude. In most cases, one of the disaggregated hazard scenarios should suffice.

Backup
==============

Because of the importance of ShakeMap, it is advisable to run redundant systems. Most 
ShakeMap operations have a primary and backup machine. The backup machine runs 
events as if it were the primary, except it does not transfer its products to the web or other 
destinations. If the primary server fails, the backup can be switched over to primary 
merely by changing the transfer configuration. This arrangement is also useful when 
software updates are available. The update can be applied and tested on the backup 
system. Once it is deemed to be operating correctly, it can be made primary, and the 
primary server can be updated.

Since most seismic networks are operated from earthquake-prone regions, there is also 
the potential that the entire facility will be taken offline. For this reason, it is desirable to 
have a backup system operating in a remote location, preferably many kilometers away.

As we have mentioned elsewhere, the USGS makes ShakeMaps for global earthquakes 
and provides backup to U.S. regional networks. If you would like to discuss remote 
backup for your ShakeMap system, please contact us.

**********************************
ShakeMap Implementation Checklist
**********************************

The checklist below is based on the one we use when discussing ShakeMap operations with active or 
potential producers within the USGS's Advanced National Seismic System (ANSS). 
While some of the issues are ANSS-specific, there may be analogous considerations for 
other regional or national networks.

1. **Triggering**

   A. Automatic Triggering System.  How is ShakeMap triggered and how does it 
      access or receive parametric data?  How is robustness of this approach 
      achieved?
   B. Location & Magnitude Reliability.  Are there limitations to location and 
      magnitude determination by the regional network that would adversely affect 
      automatic ShakeMap products? 
   C. Regional Coverage.  What are the boundaries of the area within which the 
      local network will generate ShakeMaps?
   D. Alarm Region.  For events outside ShakeMap boundaries, is a ShakeMap run 
      initiated?  Under what conditions?
   E. ShakeMap ID.  Does the naming of ShakeMap ID follow the ANSS 
      convention?  If not, can they be easily associated with the authoritative ID?

2. **Station Coverage and Parametric Data**

   A. Real-time or near--real-time data flow.  What are the types and distribution of 
      stations contributing to ShakeMap? Are all stations "ShakeMap-quality”?
   B. Parametric Data.  How are the parametric data computed? (Five parameters: 
      PGA, PGV, and three periods of PSA.) 
   C. Are parametric data imported from other sources (NSMP-triggered stations, 
      state or commercial agencies, neighboring networks, etc.)? How are these 
      integrated with the ShakeMap input?
   D. Are "Did You Feel It?" data used as input? 
   E. Co-location of different sensor types, priorities, and preventing redundant 
      input data. How are co-located instruments resolved by the network to 
      produce only a single (best) set of amplitudes for ShakeMap?

3. **System Specifications**

   A. Grind parameters. Review the parameters in *grind.conf*. How were they 
      determined?

      a. GMPEs. Which Ground-Motion Prediction Equations are used, and 
         under what conditions?
      b. IPEs. Which Intensity Prediction Equations are used, and under what 
         circumstances?
      c. GMICEs. Which Ground-Motion--Intensity Conversion Equations are 
         used?
      d. Site Amplification.  How are site conditions established and what 
         amplifications are used (GMPE-native, Borcherdt-style)?
      e. Other parameters. Grid spacing, map area, outlier levels, bias 
         parameters. Have all parameters been evaluated for optimal 
         performance?
      f. *Shake.conf*. When is map size increased, PSA and HAZUS output 
         produced, etc.?

   B. Spatial Correlation Function. Which spatial correlation function is used?
   C. Basin response. Is a basin response applied in any areas? If so, how was the 
      basin depth file produced, and are predicted ground motions consistent with 
      reality?

4. **Operations**

   A. Which version of ShakeMap is operational? Who is responsible for updating 
      the software when updates are released? When and how are the updates performed?
   B. Who is responsible for routine scientific review of ShakeMaps produced by 
      the network? Do these people receive alarms when ShakeMaps are produced?
   C. Who is responsible for routine operational review of the ShakeMap system 
      (checking logs, process and database monitoring, etc.)? When are reviews 
      performed?
   D. Reprocessing. Under what circumstances are events reprocessed (new data, 
      change in source parameters, etc.)? What about in the longer term (ShakeMap 
      software updates, changes in operational parameters)?
   E. Finite faults. For larger earthquakes, who is responsible for producing a finite 
      fault model for inclusion in ShakeMap? What procedures are in place for 
      assuring this is done?
   F. Aftershock exclusion. How will you change the triggering threshold 
      immediately after a major earthquake in your region?
   G. Version history. Under what circumstances are maps (and their input data) 
      preserved using ShakeMap versioning?
   H. Have there been any local changes to the ShakeMap software that will hinder 
      upgrades? Can these customizations be incorporated into the ShakeMap 
      distribution for easier upgrades? If not, how can they be structured to 
      accommodate easy upgrades of ShakeMap?
   I. What is the hardware for ShakeMap processing and for local web service?  
   J. How is hardware redundancy achieved?  
   K. Are the hardware and software systems automatically monitored? Do they 
      generate alerts when problems are detected?

5. **Product Distribution and Uniformity**

   A. Are products delivered to Earthquake Program Web Servers via PDL?
   B. Are local webpages produced? Where do they reside? How is ShakeMap 
      transferred? Are redundant web servers and 24/7 support available? 
   C. Are regional ShakeMap webpages customized to reflect regional 
      configurations and implementation specifics?

6. **ANSS Coordination**

   A. Provide Software/Feedback to ANSS.  To benefit current operators and to 
      ensure compatibility and ease of installing new ShakeMap software releases, 
      changes to ShakeMap software (above and beyond configuration changes) 
      should be provided to Bruce Worden for review, standardization, and 
      inclusion in new releases. 
   B. Provide contacts, their background, and roles in implementation, coordination, 
      and operations.
   C. Are all responsible parties subscribed to the *shake-dev* mailing list?

7. User Coordination:
   List significant users and outline any outreach efforts or plans. It is very useful to 
   have a feeling for which users will rely on ShakeMap in each region, as well as to 
   coordinate efforts for users of ShakeMaps for multiple regions (e.g., FEMA, 
   DHS, Military). 

8. **Scenarios and Archives**

   A. Scenario earthquakes should be made to be consistent with USGS National 
      Hazard Maps, both with attenuation relations and in source parameterization. 
      Coordination with the National Earthquake Information Center (NEIC) is essential.
   B. Is a copy of scenarios also available on the USGS web site?  
   C. How and when will scenarios be reprocessed?
   D. Archive “final" ShakeMaps for significant events.  Many users want 
      ShakeMaps for significant events "frozen in time". Once a ShakeMap gets 
      used as a reference for damage-loss modelers, insurance investigators, and 
      researchers, there needs to be an archival version of these events. Once all the 
      available ground-motion data have been collected and included in ShakeMap, 
      that Version of the map needs to be kept available even if additional updates 
      are made. (This process has not yet been fully vetted.)

9. **Backup Strategy**

   A. If the primary system fails, what provisions exist for a backup system or 
      another network to take over ShakeMap operations? Is this backup automatic 
      or manual?
   B. If the entire facility goes offline, is there an off-site backup?
   C. Are waveform or parametric data transmitted to NEIC for national-level 
      backup?

10. **Feedback**:
    Do you have any recommendations for further support, software, features, etc.? 

.. _sm35_software_guide:

*******************************************
Software Availability & Software Guide
*******************************************

ShakeMap requires the freely available PERL, MySQL, and GMT (Generic Mapping Tools), 
as well as a few other packages. PERL and GMT are used quite extensively, so any background 
with them is advantageous. You will need to assemble the basic GMT-formatted basemaps, 
road, city data files, etc., but such data may already be available for your area.

The ShakeMap software is freely available, open-source, and distributed under a Public
Domain License. It runs on Solaris, FreeBSD, Mac OS X, (U)nix, and numerous versions of Linux (including Red Hat and Debian). It 
does not run on Windows. See the Software Guide for more information. The software is available as a `SubVersion 
<http://subversion.apache.org/>`_ checkout from:

https://vault.gps.caltech.edu/repos/products/shakemap/tags/release-3.5/

.. note:: Do not attempt to install ShakeMap on Ubuntu Linux. It has been nothing 
          but a problem for everyone who has tried it, and we will no longer provide 
          support for this operating system.

The Software Guide included in the *doc* directory of the distribution will always be the 
most up-to-date and should be consulted when installing and configuring ShakeMap. The 
Software Guide may also be obtained by `download <_static/SoftwareGuideV3_5.pdf>`_.
This version of Guide is not guaranteed to be the most up-to-date, however. It should be 
used only to familiarize oneself with the general requirements of installing and operating 
ShakeMap. When installing the software, the Guide in the *doc* directory of the software
distribution should be followed.

We strongly recommend that ShakeMap operators and users sign up for the *shake-dev* mailing list:

https://geohazards.usgs.gov/mailman/listinfo/shake-dev

We use this mailing list to communicate software updates, as well as provide support 
when users have problems, suggestions, etc.

.. [#] The VS30 server currently provides GMT *grd* files in pixel node registration and 
       ShakeMap works in gridline node registration. You can fix your Vs30 file by:

       grdsample your_vs30_grid.grd -Gnew_file_name.grd –T

       You then configure *grind.conf* to look at *new_file_name.grd*. 
       See *grind.conf* for details.

