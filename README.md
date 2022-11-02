# EMTFAnalysis
A collection of scripts used in EMTF Analysis on CERN's Compact Muon Solenoid 


<h2> Introduction </h2>

These scripts use data coming from CMS's EMTF subsystem in order to measure EMTF's rates and efficiencies under a variety of conditions. 
These data come primarily from EMTFNtuples, which can be generated within the CMSSW framework, found here: https://github.com/eyigitba/EMTFTools/tree/master/EMTFNtuple.


<h3> Vocab </h3>
Let's first define some terms important to using these EMTFNtuples: 
<br><br>

**Bunch-crossing (BX)**: protons are filled in the LHC in groups known as bunches. These bunches cross CMS and collide with eachother every 25ns, and are therefore used as a time-unit known as bunch-crossings (BX). If the trigger is working correctly, muons events should be assigned BX0.

**Event**: a representation of the instantaneous state of the CMS detector at a given time. Event data is collected when a bunch-crossing is triggered upon by either CMS's muon system or ECAL system.

**emtfTrack**: anything with the prefix emtfTrack means that this parameter describes a track created by EMTF's trigger. This track is created very quickly in firmware at the time of the event, and its PT/quality are measured to determine whether the event should be triggered, and saved for analysis. These EMTF tracks are not used in final physics analysis, since they are quick and inaccurate. They only exist to determine whether a physics event is interesting (high PT usually).

**PT**: the transverse momentum of a particle, meaning the momentum of the particle in the direction perpendicular to the proton beam

**Mode**: 4-bit number representing the stations with hits in a muon track {Station 1, Station 2, Station 3, Station 4}. For example, a mode 11 (b1011) track contains hits in stations 1, 3, and 4. Tracks that have hits in at least stations 1 and 3/4 (modes 9, 10, 11, 13, 14, 15) are considered of quality 12. 

**Eta**: coordinate quanity mapped to theta, the azimuthal coordinate (considering beam as the z-axis). Used in physics contexts because of its Lorentz-Invariant(ish) properties.


<h3> EMTF Geometry </h3>

![alt text](https://twiki.cern.ch/twiki/pub/Main/NEUNewCoopGuide/MuonSystem.png)

![alt text](https://twiki.cern.ch/twiki/pub/Main/NEUNewCoopGuide/CSCLayout.png)

<h3> Efficiency: Tag and Probe </h3>

Efficiency is calculated using the tag and probe method, which is outlined as follows: 
1. Iterate through each trigger-event in the EMTFNtuple file. This file is created using muon data, meaning each event was stored because a muon was found by the trigger in EMTF.
2. Use reconstructed muons, which are more accurate at the cost of being slower to calculate. Find a reconstructed muon which should have been triggered upon (PT > 22GeV, high quality, sufficiently isolated, at BX0 (at the time of trigger)) and match it to a likewise-triggerable track that was created by EMTF. This is called the tag muon.
3. Within these events, find a second reconstructed muon called the probe, which is distinct from the tag muon. In this way we account for the biases introduced by the fact that, by definition, most events contain a high-pT triggered muon (the tag). 
4. Increment the denominator because a triggerable muon was found.
5. This probe muon, which was created in a Z->mu-mu- event, should have been triggered upon by EMTF. If it, like the tag, has a track PT > 22GeV, high-quality, and BX0 (meaning EMTF found it triggerable),then the numerator is incremented.
6. EMTF efficiency is this numerator over the denominator

We can bin this efficiency with respect to eta/phi, PT, mode, and more.

<h2> Workflow </h2>


1. Start with choosing the script you would like to run. A description is provided at the start of each script
2. Adjustable flags are also provided at the start of the script  (i.e. Do we want to require BX0? Do we want to use Ideal Geometry Data? What is the maximum number of events we want to use?). Use these flags and comments to select which input file you would like to use or add a new input path. Also define the path to your output file (make sure your input matches the forms of the alternate inputs provided).
3. **Parallelization**: If you have a lot of data, these scripts will take a long time to run. Parallelize the data in one of two ways: 
    1. `./multijob.sh <script_to_run.py> <number of jobs>` This usually works best if condor is down or if your job isn't too bulky. 20 jobs usually works well.
    2. `python condor_submit.py -i <script_to_run.py> -n <number of jobs>` Only works if you have access to condor, but much better for larger jobs. 50 jobs usually works well.
    3. This will generate <n> output files at `plots/tmp/<your_output_nameindex>.root` where index is the number identifier of your completed job. Merge these with `python recombine.py -o <your/output/file.root> -n <number of jobs>`
4. Draw the plots in CMS style to a PDF using `python <output_file_directory>/draw_plots.py` (if your output file is in one of the predefined directories). Again, make sure your output ROOT file is correctly referenced as an input in draw_plots.py. You can find the PDF's in `<output_file_directory>/pdfs/`.


