.. _sec_related-systems:

===================
Related Systems
===================
Here is a brief listing of USGS Earthquake Hazards Program rapid earthquake information products:

* `Earthquake Notification System <https://sslearthquake.usgs.gov/ens>`_ sends automated, customizable notifications of earthquakes through email, pager, or cell phone. 
* `Realtime Earthquake Map <http://earthquake.usgs.gov/earthquakes/map/>`_: Automatic maps and event information
  displayed online within minutes after earthquakes worldwidema
* `ShakeMap <http://earthquake.usgs.gov/earthquakes/shakemap/>`_ automatically generates maps displaying
  instrumentally measured shaking intensities.
* `Did You Feel It? <http://earthquake.usgs.gov/earthquakes/dyfi/>`_: Map of earthquake effects derived from citizen
  input via online Web forms
* `PAGER`_ (Prompt Assessment of Global Earthquakes for Response) rapidly
  compares the population exposed to various shaking intensities to estimate likely fatalities and economic losses
* `CISN Display <http://www.cisn.org/software/cisndisplay.html>`_: Stand-alone application that graphically alerts
  users, in near real time, of earthquakes and related hazards information.
* `ShakeCast <http://earthquake.usgs.gov/research/software/shakecast/>`_: An application for automated delivery of
  ShakeMaps and potential damage or inspection priority for specific user-selected facilities. 
* `ShakeAlert <http://www.shakealert.org/faq/>`_: Prototype Earthquake Early Warning (EEW) System.
  
While ShakeMap has met with success as a standalone product for communicating
earthquake effects to the public and the emergency response and recovery 
community, it is increasingly being incorporated into value-added products that
help in the assessment of earthquake impacts for response management and
government officials.

.. figure::  _static/SMap_SCast_DYFI_PAGER.png
   :width: 650px
   :alt: Related Systems
   :align: center
   :target: Related Systems

   Interplay between ShakeMap, DYFI, ShakeCast, and PAGER.	    

As discussed in detail the :ref:`technical-guide`, ShakeMap is augmented by
DYFI input for constraining intensities, and from those, estimates of peak 
ground motions (in some cases, and for some regions), as shown in Figure
#related-systems. DYFI and ShakeMap in conjunction then represent the shaking
hazard 
input for two other primary systems that estimated losses: ShakeCast and PAGER.
ShakeCast is intended for specific users to prioritize response for specific 
user-centric portfolios of facilities; PAGER is for more general societal-impact
assessments, providing estimated loss of life and economic impacts for the
region affected. 

.. _sec_shakecast:

ShakeCast
---------------------------------------------------

`ShakeCast`_ is a freely available
post-earthquake situational awareness application that automatically retrieves
earthquake shaking data from ShakeMap, compares intensity measures against
users’ facilities, and generates potential damage assessment notifications,
facility damage maps, and other web-based products for emergency managers and
responders.

.. sidebar:: USE CASE #2
	     
 The `California Department of Transportation
 <http://www.earthquakeauthority.com/>`_ (Caltrans) employs ShakeMap
 for post-earthquake overpass and bridge assessments for significant
 California earthquakes:

 *The Caltrans ShakeCast system performed
 reliably and as expected during the Napa earthquake. The system delivered key
 information on the potential impacts to the state bridge inventory
 within 11 minutes of the event.  Of a total of 2720 state bridges
 within the ShakeMap region, 87 state bridges were identified by
 ShakeCast as having undergone significant ground shaking and were
 assigned an inspection priority state. Of the 87 state bridges
 identified, 9 were later confirmed to have sustained minor damage.
 These 9 state bridges were ranked in the top 40% of the ShakeCast
 list.*
 --(:ref:`Turner, 2014 <turner2014>`)

ShakeCast, short for ShakeMap Broadcast, is a fully automated system for
delivering specific ShakeMap products to critical users and for triggering
established post-earthquake response protocols. ShakeMap was developed
and is used primarily for emergency response, loss estimation, and public
information; for an informed response to a serious earthquake, critical users
must go beyond just looking at ShakeMap, and understand the likely extent and
severity of impact on the facilities for which they are responsible. To this
end, the USGS has developed ShakeCast.

.. figure::  _static/Caltrans_Napa_Report.*
   :width: 350px 
   :alt: Caltrans
   :align: center
   :target: Caltrans Napa

   Caltrans ShakeCast report for the 2011 M6.0 Napa, CA earthquake. 

ShakeCast allows utilities, transportation agencies, businesses, and other
large organizations to control and optimize the earthquake information they
receive. With ShakeCast, they can automatically determine the shaking value at
their facilities, set thresholds for notification of damage states for each
facility, and then automatically notify (by pager, cell phone, or email)
specified operators and inspectors within their organizations who are
responsible for those particular facilities so they can set priorities for
response.

.. sidebar:: USE CASE #3

  *"Thought you might like to see the [Division of Safety of Dams]
  ShakeCast message for the recent Napa [Aug, 2014] Earthquake.  We have since
  divided the 1250 dams into three fragility classes (called levels of
  concern). The message provides explicit instructions on what action
  to take for each dam and transmits owner contact information. The
  message was received in my inbox 16 minutes after the earthquake,
  which was about 10 minutes after the ShakeMap version 1 was
  released. The technology has become very well accepted by the field
  inspectors. Thanks for such a great product."*
  --W. A. Fraser, C.E.G.,
  Chief, Geology Branch, CA Division of Safety of Dams, Sacramento, CA.

PAGER
---------------------------------------------------
 
Another important USGS product that uses ShakeMap output as its primary data
source is `PAGER`_ (Prompt Assessment of Global Earthquakes for Response), an
automated system that produces content concerning the impact of significant
earthquakes around the world, informing emergency responders, government and aid
agencies, and the media of the potential scope of the disaster. PAGER rapidly
assesses earthquake impacts by comparing the population exposed to each level of
shaking intensity with models of economic and fatality losses based on past
earthquakes in each country or region of the world. Earthquake alerts---which
were formerly sent based only on event magnitude and location, or population
exposure to shaking---will now be generated based also on the estimated range of
fatalities and economic losses.

PAGER alerts are based on the “Earthquake Impact Scale” developed by :ref:`Wald et al. \(2011\) <wald2011b>`.

.. figure::  _static/Nepal_M7_8_onepager_V5.*
   :width: 350px
   :alt: Nepal onePAGER 
   :align: right
   :target: Nepal OnePAGER Alert Example 

   Nepal OnePAGER Alert Example  

Public and Private Sector Tools
---------------------------------------------------
Alternatives, modifications, and enhancements to the ShakeMap methodology are
used widely around the world. Likewise, downstream derivative products and systems for loss estimation are
widely employed, both in the public and private sector. What follows is
a brief (and incomplete) description of some of these systems. Many
proprietary hazard and loss modeling systems exist in the private
sector, and typically they are openly described or referenced. 

On the shaking hazard front, domestically, some public/private sector
utilities run in-house shaking aggregation and estimation systems, 
including the East Bay Metropolitan Utility District (EBMUD's Marconi
system) and Pacific Gas and Electric (PG&E).

Impressive systems also exist in Japan, Taiwan, New Zealand, Turkey,
among other countries.

* JMA
* GNS
* INGV
  
On the rapid loss estimation front, several systems are in place in the U.S. 

Internationally, :ref:`Erdik et al. \(2011\) <erdik2011>`
and :ref:`Erdik et al. \(2014\) <erdik2014>` provide examples of
operative rapid earthquake loss estimation systems.

* Taiwan Earthquake Rapid Reporting System,
* Realtime Earthquake Assessment Disaster System in Yokohama
* Real Time Earthquake Disaster Mitigation System of the Tokyo Gas
  Co.
* IGDAS Earthquake Protection System
* Istanbul Earthquake Rapid Response System
* ELER
* SELENA
* OpenQuake (OQ, GEM Foundation)
* GDACS
* QuakeLoss (WAPMERR)
* PAGER (USGS)
  
.. note:: Links and pointers to non-USGS sites are provided for information only and do not constitute endorsement by the USGS (see `USGS policy and disclaimers <http://www.usgs.gov/laws/info_policies.html>`_).

Lastly, many systems are available and in operation in the U.S. for
aggregating hazard and impact information for emergency response and
awareness. Many are multihazard oriented, and only those with focus on
earthquake information are mentioned here. Some examples include:

* InLet (ImageCat,Inc.)
* HAZUS-MH,
* ArcGIS online

As summarized by :ref:`Gomberg and Jokobitz \(2013\) <gomberg2013>`:
“others have built in-house systems to organize, share and display observations
using commercial applications like Microsoft’s Streets and Trips and SharePoint,
Google’s GoogleEarth, or ESRI’s ArcGIS. WebEOC, a real-time Web-enabled crisis
information management system developed commercially by Esri, is meant to be an
official link among public sector emergency managers in Washington State (see
http://www.esi911.com/esi). While used by many agencies, it always was just one
of multiple communication tools. A commonly expressed desire was for a
centralized hub for all types of disaster information (like the
Department of Homeland Security’s `Virtual USA
<https://www.dropbox.com/home/Correlation/figures/SanDiego?preview=economic+losses0.png>`_)."

Further information on private sector tools can
be found in the Department of Homeland Security
(DHS) summary for the `Capstone 2014
<http://www.cusec.org/capstone14/documents/2014.03.06_PSW/2014.03.06_CAPSTONE_Private_Sector_GIS.pdf>`_
National Level (scenario) Exercise. 


