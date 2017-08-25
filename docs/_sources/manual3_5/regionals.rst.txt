.. _sec_regionals:

=======================================
 Regional Operations
=======================================


.. note::
   This section will be updated with Regional ShakeMap specifications after input from the 
   ANSS Regional Seismic Networks operators.

As described in the section on real-time ShakeMap archives, ShakeMaps are generated via independent systems running at ANSS
Regional Seismic Systems (RSNs) in Northern California, Southern California, the
Pacific Northwest, Utah, Nevada, and Alaska. For the rest of the U.S., the
ShakeMap group at the USGS National Earthquake Information Center (NEIC)
produces maps for the regional networks operating in Hawaii, New England, and
the Central and Eastern U.S. on a system referred to as Global ShakeMap (GSM).
The input, metadata, and output files produced by all these instances are
aggregated by the USGS via the Earthquake Hazards Program web system. GSM also provides
backup capabilities for the RSNs, but with degraded capabilities; not all data
are flowing from the RSNs to GSM automatically.

Separate independent systems running in Puerto Rico and New England generate
ShakeMaps, but these instances do not deliver them through the USGS Earthquake
Hazards Program webpages (at the time of this writing). GSM covers these regions, but
does not yet access the full set of data available to these regional
systems. 

In this section, we describe customizations employed by ShakeMap systems running 
throughout the Advanced National Seismic System (ANSS) regions
nationwide as well as the Global ShakeMap (GSM) system running at the NEIC in
Golden, Colorado.

.. note::
   Specifications of input parameters, data, and other configurations
   used for any ShakeMap can be found in the event-specific summary
   files (*info.xml*).

Details for about regional configurations for ShakeMap operators can
be found in the :ref:`software-guide`. For contact information for
ANSS regional operators see the :ref:`acknowledgments`. 

* **Northern California**
* **Southern California**
* **Pacific Northwest**
* **Intermountain West**
* **Mid-America**
* **Northeast**
* **Alaska**
* **Puerto Rico**
* **Hawaii**
  
  - *Coverage Area*. State of Hawaii (Bounds:)
  - *Operations*. ShakeMap operated at NEIC in conjunction with HVO RSN operators
  - *Triggering and Data Flow*.
  - *Site Condition Map*. Vs30 from topographic slope-based (:ref:`Allen and Wald, 2009b) <allen2009b>` 
  - *Ground Motion Prediction and Conversion Equations (GMPE/IPE/GMICE)*.
  - *Other Local Characteristics*.


