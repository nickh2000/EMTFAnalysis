#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src
cd /afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/ 
python high_eta_study.py -n 50 -i 3 

