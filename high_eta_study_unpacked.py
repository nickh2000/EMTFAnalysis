from cmath import tan
import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
import Helper as h
import argparse

'''
This script exists to find the efficiency of CMS's Endcap Muon Track-Finder System;

Uses a tag-and-probe method by which an event is "tagged" if a medium-quality reconstructed muon is matched to a EMTF-triggerable track.
The denominator in the efficiency ratio is the number of tagged events that also contain a distinct and medium-quality reconstructed muon
The numerator in the efficiency ratio is the number of these probe-muons that generate an EMTF-triggerable track
For most cases, we look at reconstructed probe's of a high PT (ex. 26GeV), which we should expect to be triggered upon

Efficiency in these PT regions is used to measure EMTF's ability to capture interesting physics events

***Differs from eff_EMTF.py in that this has plots to focus on the region of |eta| > 2.2, where inefficiencies were discovered

'''

#Configure input arguments
#num_jobs: split our input files into n different processes
#index: the index of the current running process, defines which chunk of input files to fetch
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--num_jobs", required=False)
parser.add_argument("-i", "--index", required = False)
args = parser.parse_args()



## Configuration settings
MAX_FILE =  -1   ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event



REQ_BX0    = True  ## Require reconstructed muon to be in BX 0  ## Require a final uGMT candidate, not just a TF muon
REQ_TIGHT = False ##Require reconstructed muon to be tight quality, not just medium
MEDIUM_ONLY = False ## Require that muons be of medium quality but not of tight quality
REQ_HLT = False ##Require muon to be included in the high-level trigger



MAX_dR  = 0.1     ## Maximum dR for L1T-offline matching
TAG_ISO = 0.15   ## Maximum relative isolation for tag muon
TAG_PT  = 26    ## Minimum offline pT for tag muon
PRB_PT  = 26    ## Minimum offline pT for probe muon
TRG_PT  = 22    ## Minimum L1T pT for probe muon
eta_min = 1.24
eta_max = 2.40
#Eta range of EMTF
modes = [range(16), [11, 13, 14, 15], [3, 5, 6, 7, 11, 13, 14, 15], [7, 11, 13, 14, 15], [7], [9], [10], [11], [13], [14], [15]] #Index plots by track modes
cuts = [22, 15, 7, 3, 0] #Index plots by probe-track pt cuts

#which modes include which station-station indices (found below)
d_phi_mode_map = {
            11:[1, 2, 5],
            13:[0, 2, 4],
            14:[0, 1, 3],
            15: [0, 1, 2, 3, 4, 5]}

#station-station indices used above
d_phi_index_map = [(1, 2), (1, 3), (1, 4),(2, 3), (2, 4), (3, 4)]

## HISTOGRAMS 
evt_tree  = TChain('EMTFNtuple/tree')

phi_mode = [[], []] #plot efficiency vs. phi for all reco muons of sufficient quality (denominator)
phi_mode_trig = [[], []] #trig plot exists to count triggered muons matched to a high-pt track for efficiency (numerator)

phi_cut = [[], []] #plot efficiency vs. phi, indexed by different trigger-pt cuts
phi_cut_trig = [[], []]


#Create, label, set range of mode-cut indexed plots
num_bins = 75
for i in range(len(modes)):
    for j in range(2):
        phi_mode[j].append(TH1D('h_phi_den_EMTF_mode_%s_%d' % (str(modes[i]), j),  '', num_bins, -3.14159, 3.14159))
        phi_mode_trig[j].append(TH1D('h_phi_num_EMTF_mode_%s_%d' % (str(modes[i]), j),  '', num_bins, -3.14159, 3.14159))
        phi_mode[j][i].GetYaxis().SetRangeUser(0, 1)
        phi_mode_trig[j][i].GetYaxis().SetRangeUser(0, 1)

#Create, label, set range of pt-cut indexed plots
for i in range (len(cuts)):
    for j in range(2):
        phi_cut[j].append(TH1D('h_phi_den_EMTF_cut_%s_%d' % (str(cuts[i]), j),  '', num_bins, -3.14159, 3.14159))
        phi_cut_trig[j].append(TH1D('h_phi_num_EMTF_cut_%s_%d' % (str(cuts[i]), j),  '', num_bins, -3.14159, 3.14159))
        phi_cut[j][i].GetYaxis().SetRangeUser(0, 1)
        phi_cut_trig[j][i].GetYaxis().SetRangeUser(0, 1)

#efficiency vs. phi in high-eta region for positive muons
h_phi_pos = TH1D('h_phi_EMTF_pos_den',  '', num_bins, -3.14159, 3.14159)
h_phi_pos_trig = TH1D('h_phi_EMTF_pos_num',  '', num_bins, -3.14159, 3.14159)

#efficiency vs. phi in high-eta region for negative muons
h_phi_neg = TH1D('h_phi_EMTF_neg_den',  '', num_bins, -3.14159, 3.14159)
h_phi_neg_trig = TH1D('h_phi_EMTF_neg_num',  '', num_bins, -3.14159, 3.14159)


if args.num_jobs:
    INDEX = int(args.index)
    NUM_JOBS = int(args.num_jobs)

if args.num_jobs:
  try: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/high_eta_study_unpacked%d.root" % (INDEX), 'create')
  except: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/high_eta_study_unpacked%d.root" % (INDEX), 'recreate')
else:
  out_file = TFile('plots/high_eta_region/high_eta_study_unpacked.root', 'recreate')
  

#Specify which datasets to use
#They must be set mutually-exclusively, if none are True, Run2 data is used
IDEAL = False
ZERO_BIAS = False
RUN2 = False
CUSTOM = False
RUN3 = True


#Determine whether to use Ideal, Run3, Run2, or customized-by-me geometry look-up tables in our EMTFNTuple Re-Emulation
#These folders are collections of EMTFNTuples, with cscSegments, re-emulated&unpacked hits, tracks enabled, ,and offline-reconstructed muons
#CMSSW Produceder for EMTFNTuples is found here: https://github.com/eyigitba/EMTFTools/tree/master/EMTFNtuple


if IDEAL:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_SingleMuon_data_13p6TeV_idealAlignment/220902_142649/0000/"
elif RUN3:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run3Alignment_2022C_v5/220925_185603/0000/"
elif CUSTOM:
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_phlutv3data_2022C_v1/221020_181346/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v9/221020_120011/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v8/221019_170559/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v7/221019_113518/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_disp_v1/221021_095909/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_init_v2/221022_102914/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_init_nonieghbor_v1/221022_103201/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_reverse_noneighbor/221023_145501/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_reverse/221023_131250/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_ref_2/221023_145705/0000/"
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_reverse_v2/221023_205023/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_ref_2_v2/221024_210356/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_ref_2_v3/221025_201619/0000/"

else:
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run2_2022C_v6/221017_155002/0000/"
    #folder = "/eos/cms/store/user/eyigitba/emtf/L1Ntuples/Run3/crabOut/SingleMuon/EMTFNtuple_Run3_SingleMuon_data_13p6TeV_355872_v1/220725_163335/0000/"
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run2Alignment_2022C_v2/220920_155151/0000/"
  
#Get all of the event files in the directory
nFiles = 0


#Flag for breaking loop if we hit max file limit
break_loop = False

#recursivelh access different subdirectories of given folder from above
for dirname, dirs, files in os.walk(folder):
    if break_loop: break

    #Get chunk from files in a given directory
if args.num_jobs and args.index:
        file_list = files[INDEX * len(files[0:MAX_FILE]) / NUM_JOBS : (INDEX + 1) * len(files[:MAX_FILE]) / NUM_JOBS]
else:
    file_list = files[0:MAX_FILE]

    #access root files in this subdirectory
for file in file_list:
    if break_loop: break
    if not '.root' in file: continue
    file_name = "%s/%s" % (dirname, file)
    nFiles   += 1
    print ('* Loading file #%s: %s' % (nFiles, file))
    evt_tree.Add(file_name)
    if nFiles == MAX_FILE: break_loop = True


#Go event by event
for event in range(evt_tree.GetEntries()):
    evt_tree.GetEntry(event)

    if event % PRT_EVT == 0:
        if args.index:
            print('high_eta_study.py: Processing Job #%d, Event #%d' % (INDEX, event))
        else: print('high_eta_study.py: Processing Event #%d' % (event))


    #List event tags and given 
    tags = []
    tracks = []


    #Need event run data to do BX cuts, because BX0 was offset to BX-1 in DAQ in Sept. 2022
    if REQ_BX0:
        run = evt_tree.eventInfo_run[0]


    #Find reconstructed (trust-worthy) muons by which to compare trigger-efficiency
    for tag in range(evt_tree.recoMuon_size):


        #Require muon tightness
        if not evt_tree.recoMuon_isMediumMuon[tag]: continue

        #Tag needs to be a certain quality to be a valid tag candidate
        if evt_tree.recoMuon_pt[tag] < TAG_PT: continue
        if evt_tree.recoMuon_iso[tag] > TAG_ISO: continue
        #Get a valid reco-muon tag coordinate
        tag_eta = evt_tree.recoMuon_etaSt2[tag]
        tag_phi = evt_tree.recoMuon_phiSt2[tag]
        if tag_eta < -99 or tag_phi < -99: 
            tag_eta = evt_tree.recoMuon_etaSt1[tag]
            tag_phi = evt_tree.recoMuon_phiSt1[tag]
        if tag_eta < -99 or tag_phi < -99: continue
        

        #Just a sanity check to clear off non-sensical muons
        if abs(tag_eta) > 2.4: continue

        #Ensure eta is valid
        if abs(tag_eta) > eta_max or abs(tag_eta) < eta_min: continue
            
        for track in range(evt_tree.emtfUnpTrack_size):
            
            track_eta = evt_tree.emtfUnpTrack_eta[track]
            track_phi = evt_tree.emtfUnpTrack_phi[track]

            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159
            
            #Ensure L1 trigger is of sufficient quality and close to the reconstructed tag
            mode = evt_tree.emtfUnpTrack_mode[track]
            if not mode in [11, 13, 14, 15]: continue

            if h.CalcDR(track_eta, track_phi, tag_eta, tag_phi ) > MAX_dR: continue
            
            if (evt_tree.emtfUnpTrack_pt[track] < TRG_PT): continue

            if REQ_BX0:
                if run >= 356798 and evt_tree.emtfUnpTrack_bx[track] != -1: continue
                elif run < 356798 and evt_tree.emtfUnpTrack_bx[track] != 0: continue

            tags.append(tag)
            tracks.append(track)
            break
        
    
    for probe in range(evt_tree.recoMuon_size):
        ####################Probe Denominator Cuts#########################################
        probe_eta = evt_tree.recoMuon_etaSt2[probe]
        probe_phi = evt_tree.recoMuon_phiSt2[probe]
        if probe_eta < -99 or probe_phi < -99: 
            probe_eta = evt_tree.recoMuon_etaSt1[probe]
            probe_phi = evt_tree.recoMuon_phiSt1[probe]
        if probe_eta < -99 or probe_phi < -99: continue 

        probe_phi_deg = probe_phi * 180. / 3.14159
        probe_pt = evt_tree.recoMuon_pt[probe]
        
        #require specific tightness
        if not evt_tree.recoMuon_isMediumMuon[probe] or abs(evt_tree.recoMuon_dxy[probe]) > .2 or abs(evt_tree.recoMuon_dz[probe]) > .2: continue

        #########Tag/Probe Isolation #####################
        matched_tag = -1
        for tag in tags:

            tag_eta = evt_tree.recoMuon_etaSt2[tag]
            tag_phi = evt_tree.recoMuon_phiSt2[tag]
            if tag_eta < -99 or tag_phi < -99: 
                tag_eta = evt_tree.recoMuon_etaSt1[tag]
                tag_phi = evt_tree.recoMuon_phiSt1[tag]

            if probe == tag: continue

            if h.CalcDR(tag_eta, tag_phi, probe_eta, probe_phi) < 4*MAX_dR: continue

            matched_tag = tag
            break
        ##################################################

        if matched_tag < 0: continue
        

        #require probe muon to be in EMTF range
        if abs(probe_eta) < eta_min or abs(probe_eta) > eta_max: continue




        #fill tightness and phi plot-denominators only with muons in the efficiency plateau, higher quality (pt > 26 GeV).ie probes that *should* be triggered upon
        if (probe_pt > PRB_PT):
            
            #limit to high eta, index by cuts
            for i in range(len(cuts)):
                if probe_eta < -2.2:
                    phi_cut[0][i].Fill(probe_phi)
                if probe_eta > 2.2:
                    phi_cut[1][i].Fill(probe_phi)

            #plot all muons of specific charge in high positive eta region
            if probe_eta > 2.2:
                if evt_tree.recoMuon_charge[probe] < 0:
                    h_phi_neg.Fill(probe_phi)
                elif evt_tree.recoMuon_charge[probe] > 0:
                    h_phi_pos.Fill(probe_phi)
            
        
    #################################################################################
        
        best_track = -1
        best_dr = -1

        ##############TRIG(Numerator) CUTS###################
        #######Match reco-probe muon to EMTF's track#############
        for track in range(evt_tree.emtfUnpTrack_size):

            track_eta = evt_tree.emtfUnpTrack_eta[track]
            track_phi = evt_tree.emtfUnpTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159
            
            new_dr = h.CalcDR( track_eta, track_phi, probe_eta, probe_phi)

            #BX Cut Trigger Track
            if REQ_BX0:
                if run >= 356798 and evt_tree.emtfUnpTrack_bx[track] != -1: continue
                elif run < 356798 and evt_tree.emtfUnpTrack_bx[track] != 0: continue
            
            #Ensure we're not matching probe muon to its tag
            if tracks[tags.index(matched_tag)] == track: continue
            
            if best_track == -1 or new_dr < best_dr: 
                best_dr = new_dr
                best_track = track


        ##############TRIG(Numerator) CUTS###################

        #make sure closest track passes distance cut
        if best_dr > 2 * MAX_dR or best_track < 0: continue
        track = best_track
        
        #Get properties of matched track
        track_eta = evt_tree.emtfUnpTrack_eta[track]
        track_phi = evt_tree.emtfUnpTrack_phi[track]
        track_phi *= 3.14159 / 180.0
        mode = evt_tree.emtfUnpTrack_mode[track]

        ##Fill denominator plots (no PT cut) for different Mode Cuts, then fill numerator plots if pt  > TRG_PT
        for i in range(len(modes)):
            if probe_pt > PRB_PT and evt_tree.emtfUnpTrack_mode[track] in modes[i]: 
                if probe_eta > 2.2:
                    phi_mode[1][i].Fill(probe_phi)
                    if evt_tree.emtfUnpTrack_pt[track] > TRG_PT:
                        phi_mode_trig[1][i].Fill(probe_phi)
                elif probe_eta < -2.2:
                    phi_mode[0][i].Fill(probe_phi)
                    if evt_tree.emtfUnpTrack_pt[track] > TRG_PT:
                        phi_mode_trig[0][i].Fill(probe_phi)
        
        ##Fill numerator plots for different PT Cuts, default high-quality mode cuts (11, 13, 14, 15)
        for i in range(len(cuts)):
            if mode in [11, 13, 14, 15] and probe_pt > PRB_PT and evt_tree.emtfUnpTrack_pt[track] > cuts[i]:
                if probe_eta < -2.2:
                    phi_cut_trig[0][i].Fill(probe_phi)
                elif probe_eta > 2.2:
                    phi_cut_trig[1][i].Fill(probe_phi)


        #Fill numerator plots for differently charge muons, standard mode and PT cuts
        if evt_tree.emtfUnpTrack_pt[track] > TRG_PT and probe_pt > PRB_PT and mode in [11, 13, 14, 15]:
            if probe_eta > 2.2:
                if evt_tree.recoMuon_charge[probe] < 0:
                    h_phi_neg_trig.Fill(probe_phi)
                elif evt_tree.recoMuon_charge[probe] > 0:
                    h_phi_pos_trig.Fill(probe_phi)

#Write Efficiency vs. phi with indexed mode cuts
for i in range(len(modes)):
    for j in range(2):
        phi_mode_trig[j][i].Write()
        phi_mode[j][i].Write()

#Write Efficiency vs. Phi indexed by trigger PT cuts
for i in range(len(cuts)):
    for j in range(2):
        phi_cut_trig[j][i].Write() 
        phi_cut[j][i].Write()

#Write efficiency plots in high-eta for negative muons
h_phi_neg_trig.Write()
h_phi_neg.Write()

#Write efficiency plots in high-eta for positive muons
h_phi_pos_trig.Write()
h_phi_pos.Write()

del out_file
