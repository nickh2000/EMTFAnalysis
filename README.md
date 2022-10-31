# EMTFAnalysis
A collection of analysis scripts used in EMTF Analysis on CERN's Compact Muon Solenoid 


<h2> Introduction </h2>

These scripts use data coming from CMS's EMTF subsystem in order to measure EMTF's rates and efficiencies under a variety of conditions. 
These data come primarily from EMTFNtuples, which can be generated within the CMSSW framework, found here: https://github.com/eyigitba/EMTFTools/tree/master/EMTFNtuple.

Let's first define some terms important to using these EMTFNtuples:

Bunch-crossing (BX): protons are filled in the LHC in groups of protons known as bunches. These bunches cross CMS every 25ns, and are therefore used as a time-unit known as bunch-crossings (BX). If the trigger is working correctly, muons events should be collected at BX0.

Event: a representation of the instantaneous state of the CMS detector at a given time. Event data is collected when a bunch-crossing is triggered upon by either CMS's muon system or ECAL system.

emtfTrack: anything with the prefix emtfTrack means that this parameter describes a track created by EMTF's trigger. This track is created very quickly in firmware at the time of the event, and its PT/quality are measured to determine whether the event should be triggered, and saved for analysis. These EMTF tracks are not used in final physics analysis, since they are quick and inaccurate. They only exist to determine whether a physics event is interesting (high PT usually).

PT: the transverse momentum of a particle, meaning the momentum of the particle in the direction perpendicular to the proton beam

Mode: 4-bit number representing the stations recording a muon hit {Station 1, Station 2, Station 3, Station 4}. For example, a mode 11 (b1011) track contains hits in stations 1, 3, and 4. Tracks that have hits in at least stations 1 and 3/4 (modes 9, 10, 11, 13, 14, 15) are considered of quality 12. 

Eta: coordinate quanity mapped to theta, the azimuthal coordinate (considering beam as the z-axis). Used in physics context's because of its Lorentz-Invariant(ish) properties.

Efficiency is calculated using the tag and probe method, which is outlined as follows: 
1. Iterate through each trigger-event in the EMTFNtuple file. This file is created using muon data, meaning each event was stored because a muon was found by the trigger in EMTF.
2. Use reconstructed muons, which are more accurate at the cost of being slower to calculate. Find a reconstructed muon which should have been triggered upon (PT > 22GeV, high quality, sufficiently isolated, at BX0 (at the time of trigger) and match it to a likewise-triggerable track that was created by EMTF. This is called the tag muon.
3. Within these events, find a second reconstructed muon called the probe, which is distinct from the tag muon. In this way we account for the biases introduced by the fact that, by definition, most events contain a high-pT triggered muon (the tag). 
4. Increment the denominator because a triggerable muon was found.
5. This probe muon, which was created in a Z->mu-mu- event, should have been triggered upon by EMTF. If it, like the tag, has a track PT > 22GeV, high-quality, and BX0 (meaning EMTF found it triggerable),then the numerator is incremented.
6. EMTF efficiency is this numerator over the denominator

We can bin this efficiency with respect to eta/phi, PT, mode, and more.
