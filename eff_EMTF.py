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

eta_min = 1.24 #Minimum eta in EMTF's range
eta_max = 2.40 #Maximum eta in EMTF's range
modes = [7, 9, 10, 11, 13, 14, 15] #Modes to be included in efficiency numerator


scale_pt_temp = [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 45, 60, 75, 100, 140, 160, 180, 200, 250, 300, 500, 1000] #high-pt range
scale_pt_temp_2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 55, 60] #low-pt range
scale_pt_2  = array('f', scale_pt_temp_2)
scale_pt  = array('f', scale_pt_temp)
max_pt = scale_pt_temp[len(scale_pt_temp) - 1] - 0.01

#pt-bounds for testing pt-resolution in different pT ranges
pt_ranges = [(0, 5), (5, 10), (10, 20), (20, 50), (0, 50)]
eta_ranges = [(1.2, 1.6), (1.6, 2.1,), (2.1, 2.4)]


## HISTOGRAMS 
evt_tree  = TChain('EMTFNtuple/tree')


#Plot efficiency vs. PT
#_trg plots represent the efficiency numerator (probe muons that were triggered upon)
h_pt  = TH1D('h_pt_den_EMTF',  '', len(scale_pt_temp_2) - 1,  scale_pt_2)
h_pt_trg = TH1D('h_pt_num_EMTF',  '', len(scale_pt_temp_2)-1,  scale_pt_2)

#Include higher PT's
h_pt_large = TH1D('h_pt_den_EMTF_large',  '', len(scale_pt_temp) - 1,  scale_pt)
h_pt_trg_large = TH1D('h_pt_num_EMTF_large',  '', len(scale_pt_temp) - 1,  scale_pt)


#Plot efficiency vs. eta
h_eta  = TH1D('h_eta_den_EMTF',  '', 50,  -2.5, 2.5)
h_eta_trg = TH1D('h_eta_num_EMTF',  '', 50, -2.5, 2.5)


#Plot color map of efficiecy at pt-eta
h_eta_pt = TH2D('h_eta_pt_den_EMTF', '', 50, -2.4, 2.4, 50, 0, 100)
h_eta_phi = TH2D('h_eta_phi_den_EMTF', '', 50, -2.4, 2.4, 30, -3.14159, 3.14159)
h_eta_pt_trg = TH2D('h_eta_pt_num_EMTF', '', 50, -2.4, 2.4, 50, 0, 100)
h_eta_phi_trg = TH2D('h_eta_phi_num_EMTF', '', 50, -2.4, 2.4, 30, -3.14159, 3.14159)

#Do same but only with tight tags
h_eta_isTight =TH2D('h_eta_tight_den_EMTF', '', 50, -2.4, 2.4, 2, 0, 2)
h_eta_isTight_trg =TH2D('h_eta_tight_num_EMTF', '', 50, -2.4, 2.4, 2, 0, 2)
h_eta_isTight_denom =TH2D('h_eta_tight_denom_EMTF', '', 50, -2.4, 2.4, 2, 0, 2) #plot percent of tight / medium muons by eta, need to divide by totoal muons of each tightness

pos_end_mode = TH2D('pos_end_mode', '', 30, -3.14159, 3.14159, 16, 0, 16) ##plot region in eta > 2.3 with low efficiency

h_eta_qual = TH2D('h_eta_qual_EMTF', '', 50, -2.4, 2.4, 16, 0, 16) #plot the mode vs eta


#Total muons of tightness and not-tightness
tight_0_tot =0 
tight_1_tot = 0


#split up input-files to parallelize code in condor_submit.py or ./multi-job.bash
if args.num_jobs:
    INDEX = int(args.index)
    NUM_JOBS = int(args.num_jobs)

#if we are parallelizing, use temporary chunked-outfiles
if args.num_jobs:
  try: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/eff_EMTF_custom_chamberTwo%d.root" % (INDEX), 'recreate')
  except: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/eff_EMTF_custom_chamberTwo%d.root" % (INDEX), 'recreate')
else:
  out_file = TFile('plots/alignment/alignment_study_single.root', 'recreate')


#Specify which datasets to use
#They must be set mutually-exclusively, if none are True, Run2 data is used
IDEAL = False
ZERO_BIAS = False
RUN2 = False
CUSTOM = True
RUN3 = False


#Determine whether to use Ideal, Run3, Run2, or customized-by-me geometry look-up tables in our EMTFNTuple Re-Emulation
#These folders are collections of EMTFNTuples, with cscSegments, re-emulated&unpacked hits, tracks enabled, ,and offline-reconstructed muons
#CMSSW Produceder for EMTFNTuples is found here: https://github.com/eyigitba/EMTFTools/tree/master/EMTFNtuple


if IDEAL: 
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_SingleMuon_data_13p6TeV_idealAlignment/220902_142649/0000/"
elif RUN3:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run3Alignment_2022C_v5/220925_185603/0000/"
elif RUN2:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run2Alignment_2022C_v5/221011_115700/0000/"
elif CUSTOM:
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v5/221017_105731/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_reverse_v4/221026_174204/0000/"
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v10/221027_083620/0000/"
    
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v13/221031_221428/0000/"
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_chamber_v3/221109_220323/0000/"

print(folder)
#Get all of the event files in the directory
nFiles = 0

files = Popen(['ls', folder], stdout=PIPE).communicate()[0].split()

#If using parallelization, get chunk of input files to use for this batch
if args.num_jobs and args.index:
    file_list = files[INDEX * len(files[:MAX_FILE]) / NUM_JOBS : (INDEX + 1) * len(files[:MAX_FILE]) / NUM_JOBS]
else:
    file_list = files[0:MAX_FILE]

#Add input files to event-tree
for file in file_list:
    if not '.root' in file: continue
    file_name = "%s%s" % (folder, file)
    nFiles   += 1
    print ('* Loading file #%s: %s' % (nFiles, file))
    evt_tree.Add(file_name)
    if nFiles == MAX_FILE: break

#Go event by event
for event in range(evt_tree.GetEntries()):
    evt_tree.GetEntry(event)

    if event % PRT_EVT == 0:
        if args.index:
            print('eff_EMTF.py: Processing Job #%d, Event #%d' % (INDEX, event))
        else: print('eff_EMTF.py: Processing Event #%d' % (event))
    

    #List event tags and given 
    tags = []
    tracks = []


    #Need event run data to do BX cuts, because BX0 was offset to BX-1 in DAQ in Sept. 2022
    if REQ_BX0:
        run = evt_tree.eventInfo_run[0]


    #Find reconstructed (trust-worthy) muons by which to compare trigger-efficiency
    for tag in range(evt_tree.recoMuon_size):
        

        #Require tag-muon tightness
        if not REQ_TIGHT and not evt_tree.recoMuon_isMediumMuon[tag]: continue
        if REQ_TIGHT and not evt_tree.recoMuon_isTightMuon[tag]: continue


        #Tag needs to be a certain quality to be a valid tag candidate
        if evt_tree.recoMuon_pt[tag] < TAG_PT: continue
        if evt_tree.recoMuon_iso[tag] > TAG_ISO: continue
        if REQ_HLT and evt_tree.recoMuon_hlt_isomu[tag] != 1:   continue
        if REQ_HLT and evt_tree.recoMuon_hlt_isoDeltaR[tag] > TAG_ISO: continue

        #Reco-muons are assigned -99 coordinates if a hit is not found at that station
        tag_eta = evt_tree.recoMuon_etaSt2[tag]
        tag_phi = evt_tree.recoMuon_phiSt2[tag]
        if tag_eta < -99 or tag_phi < -99: 
            tag_eta = evt_tree.recoMuon_etaSt1[tag]
            tag_phi = evt_tree.recoMuon_phiSt1[tag]
        if tag_eta < -99 or tag_phi < -99: continue


        #Just a sanity check to clear off non-sensical muons
        if abs(tag_eta) > 2.4: continue


        #We need to match the tag to the track by which EMTF found it (this is a pre-biased Muon data-set)
        for track in range(evt_tree.emtfTrack_size):
            
            track_eta = evt_tree.emtfTrack_eta[track]
            track_phi = evt_tree.emtfTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159


            if REQ_BX0:
                if run >= 356798 and evt_tree.emtfTrack_bx[track] != -1: continue
                elif run < 356798 and evt_tree.emtfTrack_bx[track] != 0: continue
            
            #Ensure L1 trigger is of sufficient quality and close to the reconstructed tag
            if (evt_tree.emtfTrack_mode[track] < 11 or evt_tree.emtfTrack_mode[track] == 12): continue
            if h.CalcDR( track_eta, track_phi, tag_eta, tag_phi ) > MAX_dR: continue
            if (evt_tree.emtfTrack_pt[track] < TRG_PT): continue

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
        probe_pt = min(probe_pt, max_pt)
        
        #require specific tightness
        if not REQ_TIGHT and not evt_tree.recoMuon_isMediumMuon[probe]: continue
        if REQ_TIGHT and not evt_tree.recoMuon_isTightMuon[probe]: continue
        if MEDIUM_ONLY and evt_tree.recoMuon_isTightMuon[probe]: continue

        if abs(probe_eta) > 2.4: continue

        #########Tag/Probe Isolation #####################
        matched_tag = -1
        for tag in tags:

            tag_eta = evt_tree.recoMuon_etaSt2[tag]
            tag_phi = evt_tree.recoMuon_phiSt2[tag]
            if tag_eta < -99 or tag_phi < -99: 
                tag_eta = evt_tree.recoMuon_etaSt1[tag]
                tag_phi = evt_tree.recoMuon_phiSt1[tag]

            if probe == tag: continue
            #Determine whether to use Ideal, Run3, Run2, or customized-by-me geometry look-up tables in our EMTFNTuple Re-Emulation
            if h.CalcDR(tag_eta, tag_phi, probe_eta, probe_phi) < 4*MAX_dR: continue

            matched_tag = tag
            break
        ##################################################

        if matched_tag < 0: continue
        

        #require probe muon to be in EMTF range
        if abs(probe_eta) < eta_min or abs(probe_eta) > eta_max: continue


        #fill pt plots with all denominators
        h_pt.Fill(probe_pt)
        h_pt_large.Fill(probe_pt)
        h_eta_pt.Fill(probe_eta, probe_pt)


        #fill tightness and phi plots only with muons in the efficiency plateau, higher quality (pt > 26 GeV).ie probes that *should* be triggered upon
        if (probe_pt > PRB_PT):

            h_eta.Fill(probe_eta)
            h_eta_isTight.Fill(probe_eta, int(evt_tree.recoMuon_isTightMuon[probe]))
            h_eta_phi.Fill(probe_eta, probe_phi)
            
        

    #################################################################################
        
        best_track = -1
        best_dr = -1

    ##############TRIG(Numerator) CUTS###################
        #######Match reco-probe muon to EMTF's track#############
        for track in range(evt_tree.emtfTrack_size):

            track_eta = evt_tree.emtfTrack_eta[track]
            track_phi = evt_tree.emtfTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159

            new_dr = h.CalcDR( track_eta, track_phi, probe_eta, probe_phi)

            #BX Cut Trigger Track
            if REQ_BX0:
                if run >= 356798 and evt_tree.emtfTrack_bx[track] != -1: continue
                elif run < 356798 and evt_tree.emtfTrack_bx[track] != 0: continue
            
            #Ensure we're not matching probe muon to its tag
            if tracks[tags.index(matched_tag)] == track: continue

            if best_track == -1 or new_dr < best_dr:
                best_dr = new_dr
                best_track = track

        #make sure closest track passes distance cut
        if best_dr > 2 * MAX_dR or best_track < 0: continue
        track = best_track


        #Plot efficiency in the high-eta region (further explored in high_eta_study.py)
        if probe_eta >= 2.3 and evt_tree.emtfTrack_pt[track] > TRG_PT: pos_end_mode.Fill(probe_phi, evt_tree.emtfTrack_mode[track]) #find incidence of trigger muons in erroneous high eta region


        #plot efficiency in space
 
        #Tag and probe momentum/quality cuts for all subsequenct plots
        if (evt_tree.emtfTrack_mode[track] < 11 or evt_tree.emtfTrack_mode[track] == 12): continue #standard single-muon quality cuts
        if evt_tree.emtfTrack_pt[track] < TRG_PT - .01: continue #single L1 trigger pt cut
        ###########################################
        
        #Fill numerator with all muons that were triggered
        h_pt_trg.Fill(probe_pt)
        h_pt_trg_large.Fill(probe_pt)
        h_eta_pt_trg.Fill(probe_eta, probe_pt)

        ##Only plot muons who are in the efficiency plateau (pt > 26 GeV)
        if (probe_pt > PRB_PT):

            #count total number of tight / loose muons
            if evt_tree.recoMuon_isTightMuon[probe]: tight_1_tot += 1
            else: tight_0_tot += 1

            #plot tightness of muons in the plateau region

            h_eta_trg.Fill(probe_eta)
            h_eta_isTight_trg.Fill(probe_eta, evt_tree.recoMuon_isTightMuon[probe])
            h_eta_phi_trg.Fill(probe_eta, probe_phi) #looser pt constraints (7GeV) and loser mode constraints to get more detail on phi plot
            d_phi_12 = list(evt_tree.emtfTrack_ptLUT_deltaPh)[track][0]

        h_eta_qual.Fill(probe_eta, evt_tree.emtfTrack_mode[track]) #fill mode plot without regard to mode cuts

print ('%s: %.1f +/- %.1f%%' % ( "EMTF", 100 * h_pt_trg.Integral() / h_pt.Integral(), 
                                    (100 * h_pt_trg.Integral() / h_pt.Integral()) * 
                                    math.sqrt(h_pt.Integral()) / h_pt.Integral() ))

#Statistical error graph of pt efficiencies


h_pt_trg.Write() 
h_pt_trg_large.Write()

h_pt_trg.Divide(h_pt)
h_pt_trg_large.Divide(h_pt_large)

h_pt_trg.SetName(h_pt_trg.GetName()+'_eff')
h_pt_trg_large.SetName(h_pt_trg_large.GetName() + "_eff")

h_pt.Write()
h_pt_large.Write()


#2D plot is similar to above, includes eta axis
h_eta_pt_trg.Write() 
h_eta_pt_trg.Divide(h_eta_pt)
h_eta_pt_trg.SetName(h_eta_pt_trg.GetName()+'_eff')
h_eta_pt_trg.Write()
h_eta_pt.Write()


#Plot efficiency of eta with respect to phi (only in 26GeV region)
h_eta_phi_trg.Write() 
h_eta_phi_trg.Divide(h_eta_phi)
h_eta_phi_trg.SetName(h_eta_phi_trg.GetName()+'_eff')

h_eta_phi_trg.Draw("colz")

h_eta_phi_trg.SetLineWidth(1)
h_eta_phi_trg.SetLineColor(1)

h_eta_phi_trg.Write()
h_eta_phi.Write()

h_eta.Write()
h_eta_trg.Write() 
h_eta_trg.Divide(h_eta)
h_eta_trg.SetName(h_eta_trg.GetName()+'_eff')
h_eta_trg.Write()
h_eta.Write()

#Save number of tight muuons
h_eta_isTight_trg.Write()


for i in range(1, 51):
    h_eta_isTight_denom.SetBinContent(i, 1, tight_0_tot)
    h_eta_isTight_denom.SetBinContent(i, 2, tight_1_tot)

h_eta_isTight_num = h_eta_isTight_trg.Clone("h_eta_tight_num") #clone triggered muon count and divide by total number of triggered muons (per tightness)
h_eta_isTight_num.Divide(h_eta_isTight_denom)
h_eta_isTight_denom.Write()
h_eta_isTight_num.Write()



#plot the muon efficiencies vs tightness
h_eta_isTight_trg.Divide(h_eta_isTight)
h_eta_isTight_trg.SetName(h_eta_isTight_trg.GetName()+'_eff')
h_eta_isTight_trg.Write()
h_eta_isTight.Write()
#plot mode vs eta
h_eta_qual.Write()

#plot mode vs phi in the high positive eta region
pos_end_mode.Write()

del out_file