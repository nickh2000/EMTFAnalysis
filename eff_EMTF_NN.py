import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
import Helper as h


## Configuration settings
USE_EMUL = True    ## Use emulated L1T muons instead of unpacked
MAX_FILE = -1       ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event

REQ_BX0    = True  ## Require L1T muon to be in BX 0  ## Require a final uGMT candidate, not just a TF muon
## REQ_Z      = False ## Require tag and probe muon to satisfy 81 < mass < 101 GeV (not yet implemented - AWB 25.04.2019)
REQ_TIGHT = False
REQ_HLT = False

MAX_dR  = 0.1     ## Maximum dR for L1T-offline matching
TAG_ISO = 0.15   ## Maximum relative isolation for tag muon
TAG_PT  = 26    ## Minimum offline pT for tag muon
PRB_PT  = 26    ## Minimum offline pT for probe muon
TRG_PT  = 22    ## Minimum L1T pT for probe muon
MIN_MODE = 12    ## Minimum L1T qual for probe muon

eta_min = 1.24
eta_max = 2.40
max_pt = 50


scale_pt_temp = [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 45, 60, 75, 100, 140, 160, 180, 200, 250, 300, 500, 1000]
scale_pt_temp_2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 55, 60]
scale_pt_2  = array('f', scale_pt_temp_2)
scale_pt  = array('f', scale_pt_temp)
max_pt = scale_pt_temp[len(scale_pt_temp) - 1] - 0.01


## L1NTuple branches
evt_tree  = TChain('EMTFNtuple/tree')
h_pt  = TH1D('h_pt_EMTF',  '', len(scale_pt_temp_2)-1,  scale_pt_2)
h_pt_trg = TH1D('h_pt_trg_EMTF',  '', len(scale_pt_temp_2)-1,  scale_pt_2)

h_dxy  = TH1D('h_dxy_EMTF',  '', 30, 0, 50)
h_dxy_trg = TH1D('h_dxy_trg_EMTF',  '', 30, 0, 50)

folder = "/eos/cms/store/user/eyigitba/emtf/L1Ntuples/Run3/crabOut/SingleMuon/EMTFNtuple_Run3_SingleMuon_data_13p6TeV_355872_v1/220725_163335/0000/"
out_file = TFile('plots/L1T_eff_Pt_EMTF_NN.root', 'recreate')

#Get all of the event files in the directory
nFiles = 0
for file in Popen(['ls', folder], stdout=PIPE).communicate()[0].split():
    if not '.root' in file: continue
    file_name = "%s%s" % (folder, file)
    nFiles   += 1
    print ('* Loading file #%s: %s' % (nFiles, file))
    evt_tree.Add(file_name)
    if nFiles == MAX_FILE: break

for event in range(evt_tree.GetEntries()):
    evt_tree.GetEntry(event)

    if event % PRT_EVT == 0:
        print('Processing Event #%d' % (event))

    tags = []
    tracks = []
    for tag in range(evt_tree.recoMuon_size):

        if not REQ_TIGHT and not evt_tree.recoMuon_isMediumMuon[tag]: continue
        if REQ_TIGHT and not evt_tree.recoMuon_isTightMuon[tag]: continue

        if evt_tree.recoMuon_pt[tag] < TAG_PT: continue
        if evt_tree.recoMuon_iso[tag] > TAG_ISO: continue
        if REQ_HLT and evt_tree.recoMuon_hlt_isomu[tag] != 1:   continue
        if REQ_HLT and evt_tree.recoMuon_hlt_isoDeltaR[tag] > TAG_ISO: continue

        tag_eta = evt_tree.recoMuon_etaSt2[tag]
        tag_phi = evt_tree.recoMuon_phiSt2[tag]
        if tag_eta < -99 or tag_phi < -99: 
            tag_eta = evt_tree.recoMuon_etaSt1[tag]
            tag_phi = evt_tree.recoMuon_phiSt1[tag]
        if tag_eta < -99 or tag_phi < -99: continue

        if abs(tag_eta) > 2.4: continue

        for track in range(evt_tree.emtfTrack_size):
            
            
            track_eta = evt_tree.emtfTrack_eta[track]
            track_phi = evt_tree.emtfTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159

            if (evt_tree.emtfTrack_mode[track] < 11 or evt_tree.emtfTrack_mode[track] == 12): continue
            if h.CalcDR( track_eta, track_phi, tag_eta, tag_phi ) > MAX_dR: continue
            if (evt_tree.emtfTrack_pt_dxy[track] < TRG_PT): continue

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

        probe_pt = evt_tree.recoMuon_pt[probe]
        probe_pt = min(probe_pt, max_pt)
        
        if not REQ_TIGHT and not evt_tree.recoMuon_isMediumMuon[probe]: continue
        if REQ_TIGHT and not evt_tree.recoMuon_isTightMuon[probe]: continue
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

            if h.CalcDR(tag_eta, tag_phi, probe_eta, probe_phi) < 4*MAX_dR: continue

            matched_tag = tag
            break
        ##################################################

        if matched_tag < 0: continue
            
        if abs(probe_eta) < eta_min or abs(probe_eta) > eta_max: continue
        h_pt.Fill(probe_pt)
        if probe_pt > PRB_PT: h_dxy.Fill(evt_tree.recoMuon_dxy[probe] * 100)

    #################################################################################
        

        for track in range(evt_tree.emtfTrack_size):

            track_eta = evt_tree.emtfTrack_eta[track]
            track_phi = evt_tree.emtfTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159


            ##############TRIG(Numerator) CUTS###################
            if tracks[tags.index(matched_tag)] == track: continue

            if h.CalcDR( track_eta, track_phi, probe_eta, probe_phi ) > 2*MAX_dR: continue

            if evt_tree.emtfTrack_pt_dxy[track] < TRG_PT - .01: continue

            if (evt_tree.emtfTrack_mode[track] < 11 or evt_tree.emtfTrack_mode[track] == 12): continue

            if REQ_BX0 and evt_tree.emtfTrack_bx[track] != 0: continue
            ###########################################
            h_pt_trg.Fill(probe_pt)
            if probe_pt > PRB_PT: h_dxy_trg.Fill(evt_tree.recoMuon_dxy[probe] * 100)
            break



print ('%s: %.1f +/- %.1f%%' % ( "EMTF", 100 * h_pt_trg.Integral() / h_pt.Integral(), 
                                    (100 * h_pt_trg.Integral() / h_pt.Integral()) * 
                                    math.sqrt(h_pt.Integral()) / h_pt.Integral() ))

h_err_pt = TGraphAsymmErrors(h_pt_trg, h_pt)
h_err_pt.SetName('h_pt_trg_EMTF_err',)
h_err_pt.Write()

h_pt_trg.Write() 
h_pt_trg.Divide(h_pt)
h_pt_trg.SetName(h_pt_trg.GetName()+'_eff')
h_pt_trg.Write()

h_err_dxy = TGraphAsymmErrors(h_dxy_trg, h_dxy)
h_err_dxy.SetName('h_dxy_trg_EMTF_err',)
h_err_dxy.Write()

h_dxy_trg.Write() 
h_dxy_trg.Divide(h_dxy)
h_dxy_trg.SetName(h_dxy_trg.GetName()+'_eff')
h_dxy_trg.Write()
del out_file