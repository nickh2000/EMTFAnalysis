#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src
cd /afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/ 
python station_station_phi.py -n 50 -i 13 

