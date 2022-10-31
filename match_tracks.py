import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
import Helper as h

## Configuration settings
MAX_FILE =  -1       ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event

'''This script investigates eta mis-matches between re-emulated and firmware tracks in EMTF'''

#Set up NTuple object and i/o files
evt_tree  = TChain('EMTFNtuple/tree')
folder = "/eos/cms/store/user/eyigitba/emtf/L1Ntuples/Run3/EfficiencyStudies/357479/EMTFNtuples/Muon/EMTFNtuple_Muon_run3BDT_357479_Run2Alignment_20220825_114518/220825_164535/0000/"
out_file = TFile('plots/L1T_eff_Pt_EMTF.root', 'recreate')

#Get all of the event files in the directory
nFiles = 0
for file in Popen(['ls', folder], stdout=PIPE).communicate()[0].split():
    if not '.root' in file: continue
    file_name = "%s%s" % (folder, file)
    nFiles   += 1
    print ('* Loading file #%s: %s' % (nFiles, file))
    evt_tree.Add(file_name)
    if nFiles == MAX_FILE: break


#Loop through all events in event-tree
mismatched_tracks = 0
num_tracks = 0
for event in range(evt_tree.GetEntries()):
    evt_tree.GetEntry(event)
    if event % PRT_EVT == 0:
        print('Processing Event #%d' % (event))
    #Count one track per event
    num_tracks += evt_tree.emtfTrack_size
    for i in range(evt_tree.emtfTrack_size):
        for j in range(evt_tree.emtfUnpTrack_size):
            #Match unpacked and re-emulated tracks by mode, phi, sector, and BX
            if evt_tree.emtfTrack_mode[i] == evt_tree.emtfUnpTrack_mode[j]:
                if evt_tree.emtfTrack_phi[i] == evt_tree.emtfUnpTrack_phi[j]:
                    if evt_tree.emtfTrack_sector[i] == evt_tree.emtfUnpTrack_sector[j]:
                        if evt_tree.emtfTrack_bx[i] == evt_tree.emtfUnpTrack_bx[j]:

                            #Consider tracks mismatched if eta differs and print its information
                            if abs(evt_tree.emtfTrack_eta[i] - evt_tree.emtfUnpTrack_eta[j]) >= .001:
                                mismatched_tracks += 1
                                print("Track 1: Mode=%f Eta=%f Phi=%f Pt=%f" % (evt_tree.emtfTrack_mode[i], evt_tree.emtfTrack_eta[i], evt_tree.emtfTrack_phi[i], evt_tree.emtfTrack_pt[i]))
                                print("Track 2: Mode=%f Eta=%f Phi=%f Pt=%f" % (evt_tree.emtfUnpTrack_mode[j], evt_tree.emtfUnpTrack_eta[j], evt_tree.emtfUnpTrack_phi[j], evt_tree.emtfUnpTrack_pt[j]))

#Get estimate of proportion of tracks that are mismatched
print("Total Tracks: " + str(num_tracks) + ", Mismatched Tracks: " + str(mismatched_tracks))
