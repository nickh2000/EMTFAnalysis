#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src
cd /afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/ 
python 2018compare.py -n 50 -i 25 

