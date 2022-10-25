from cmath import tan
import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
import Helper as h
import argparse

def cms_latex():
  cms_label = TLatex()
  cms_label.SetTextSize(0.04)
  cms_label.DrawLatexNDC(0.1, 0.92, "#bf{ #font[22]{CMS} #font[72]{Efficiency Studies}}");
  return cms_label


def head():
  header = TLatex()
  header.SetTextSize(0.03)
  header.DrawLatexNDC(0.63, 0.92, "#sqrt{s} = 13.6 TeV, Run 3 Data");
  return header

'''This script exists to analyze tagged events in the high-eta regions of EMTF's positive endcap'''

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--num_jobs", required=False)
parser.add_argument("-i", "--index", required = False)
args = parser.parse_args()



## Configuration settings
MAX_FILE =  -1     ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event

IDEAL = False
RUN3 = False
CUSTOM = True

REQ_BX0    = True  ## Require L1T muon to be in BX 0  ## Require a final uGMT candidate, not just a TF muon
REQ_TIGHT = True 
MEDIUM_ONLY = False ## Require that muons be of medium quality but not of tight quality


MAX_dR  = 0.1     ## Maximum dR for L1T-offline matching
TAG_ISO = 0.15   ## Maximum relative isolation for tag muon
TAG_PT  = 26    ## Minimum offline pT for tag muon
PRB_PT  = 26    ## Minimum offline pT for probe muon
TRG_PT  = 22    ## Minimum L1T pT for probe muon

#Eta range of EMTF
eta_min = 1.24
eta_max = 2.40
max_pt = 50
modes = [range(16), [11, 13, 14, 15], [7, 9, 10, 11, 13, 14, 15], [7], [9], [10], [11], [13], [14], [15]] #Index plots by track modes
cuts = [22, 15, 7, 3] #Index plots by probe-track pt cuts

#plot specific pt's
scale_pt_temp = [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 45, 60, 75, 100, 140, 160, 180, 200, 250, 300, 500, 1000]
scale_pt_temp_2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 55, 60]
scale_pt_2  = array('f', scale_pt_temp_2)
scale_pt  = array('f', scale_pt_temp)
max_pt = scale_pt_temp[len(scale_pt_temp) - 1] - 0.01

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


#plot a track's delta-phi between stations as a function of phi
#encap, 2-station-tuple, sector (0 index is all sectors)
d_phi_plots = np.zeros(shape = (2, 6, 7, 2), dtype=object)

num_bins = 50
for i in range(len(modes)):
    for j in range(2):
        phi_mode[j].append(TH1D('h_phi_den_EMTF_mode_%s_%d' % (str(modes[i]), j),  '', num_bins, -3.14159, 3.14159))
        phi_mode_trig[j].append(TH1D('h_phi_num_EMTF_mode_%s_%d' % (str(modes[i]), j),  '', num_bins, -3.14159, 3.14159))
        phi_mode[j][i].GetYaxis().SetRangeUser(0, 1)
        phi_mode_trig[j][i].GetYaxis().SetRangeUser(0, 1)

for i in range (len(cuts)):
    for j in range(2):
        phi_cut[j].append(TH1D('h_phi_den_EMTF_cut_%s_%d' % (str(cuts[i]), j),  '', num_bins, -3.14159, 3.14159))
        phi_cut_trig[j].append(TH1D('h_phi_num_EMTF_cut_%s_%d' % (str(cuts[i]), j),  '', num_bins, -3.14159, 3.14159))
        phi_cut[j][i].GetYaxis().SetRangeUser(0, 1)
        phi_cut_trig[j][i].GetYaxis().SetRangeUser(0, 1)


for endcap in range(2):
    for d in range(6):
        for charge in range(2):
            #index 0 includes all sectors to get a global dphi offset
            d_phi_plots[endcap][d][0][charge] = TH1D('d_phi_EMTF_%s_%d_c%d_all' % (str(d_phi_index_map[d]), endcap, charge),  '', 128, -64, 64)

            #then plot by sector
            for sector in range(6):
                d_phi_plots[endcap][d][sector + 1][charge] = TH1D('d_phi_EMTF_%s_%d_%d_c%d' % (str(d_phi_index_map[d]), endcap, sector + 1, charge),  '', 128, -64, 64)
                d_phi_plots[endcap][d][sector + 1][charge].GetYaxis().SetTitle('Instances')
                d_phi_plots[endcap][d][sector + 1][charge].GetXaxis().SetTitle("#Delta#phi_{%d%d}" % (d_phi_index_map[d][0], d_phi_index_map[d][1]))

                d_phi_plots[endcap][d][sector + 1][charge].SetTitle("#Delta#phi_{%d%d}, Endcap %d, Sector %d" %
                (d_phi_index_map[d][0], 
                d_phi_index_map[d][1],
                1 if endcap else -1, 
                sector + 1))

#efficiency vs. phi for positive muons
h_phi_pos = TH1D('h_phi_EMTF_pos_den',  '', num_bins, -3.14159, 3.14159)
h_phi_pos_trig = TH1D('h_phi_EMTF_pos_num',  '', num_bins, -3.14159, 3.14159)

#efficiency vs. phi for negative muons
h_phi_neg = TH1D('h_phi_EMTF_neg_den',  '', num_bins, -3.14159, 3.14159)
h_phi_neg_trig = TH1D('h_phi_EMTF_neg_num',  '', num_bins, -3.14159, 3.14159)

#encap, station, ring
phi_diff_plot = np.zeros(shape=(2, 4, 3), dtype=object)

for endcap in range(2):
    for station in range(4):
        for ring in range(3):
                phi_diff_plot[endcap][station][ring] = TH2D('phi_diff_EMTF_%d_%d_%d' % (endcap, station + 1, ring + 1), '', 36, 1, 37, 100, -.5, .5)
                
                

if args.index:
    INDEX = int(args.index)
    NUM_JOBS = int(args.num_jobs)

if args.num_jobs:
  try: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/high_eta_study_unpacked%d.root" % (INDEX), 'create')
  except: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/high_eta_study_unpacked%d.root" % (INDEX), 'recreate')
else:
  out_file = TFile('plots/high_eta_region/high_eta_study_unpacked.root', 'recreate')


if IDEAL:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_SingleMuon_data_13p6TeV_idealAlignment/220902_142649/0000/"
elif RUN3:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run3Alignment_2022C_v5/220925_185603/0000/"
elif CUSTOM:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v5/221017_105731/0000/"
else:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run2Alignment_2022C_v2/220920_155151/0000/"
  

#Get all of the event files in the directory
nFiles = 0

files = Popen(['ls', folder], stdout=PIPE).communicate()[0].split()

if args.num_jobs and args.index:
      file_list = files[INDEX * len(files[:MAX_FILE]) / NUM_JOBS : (INDEX + 1) * len(files[:MAX_FILE]) / NUM_JOBS]
else:
    file_list = files[0:MAX_FILE]

for file in file_list:
    if not '.root' in file: continue
    file_name = "%s%s" % (folder, file)
    nFiles   += 1
    print ('* Loading file #%s: %s' % (nFiles, file))
    evt_tree.Add(file_name)

for event in range(evt_tree.GetEntries()):
    evt_tree.GetEntry(event)
    if event == MAX_EVT: break
    if event % PRT_EVT == 0:
        if args.index:
            print('high_eta_study_unpacked.py: Processing Job #%d, Event #%d' % (INDEX, event))
        else: print('high_eta_study_unpacked.py: Processing Event #%d' % (event))

    tags = []
    tracks = []

    run = evt_tree.eventInfo_run[0]
    for tag in range(evt_tree.recoMuon_size):


        #Require muon tightness
        if not REQ_TIGHT and not evt_tree.recoMuon_isMediumMuon[tag]: continue
        if REQ_TIGHT and not evt_tree.recoMuon_isTightMuon[tag]: continue


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
        

        #Ensure eta is valid
        if abs(tag_eta) > 2.4: continue

        for track in range(evt_tree.emtfUnpTrack_size):
            
            
            track_eta = evt_tree.emtfUnpTrack_eta[track]
            track_phi = evt_tree.emtfUnpTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159
            
            #Ensure L1 trigger is of sufficient quality and close to the reconstructed tag
            if (evt_tree.emtfUnpTrack_mode[track] < 11 or evt_tree.emtfUnpTrack_mode[track] == 12): continue
            if h.CalcDR(track_eta, track_phi, tag_eta, tag_phi ) > MAX_dR: continue
            if (evt_tree.emtfUnpTrack_pt[track] < TRG_PT): continue
            if run >= 356798 and evt_tree.emtfUnpTrack_bx[track] != -1: continue
            elif run < 356798 and evt_tree.emtfUnpTrack_bx[track] != 0: continue

            tags.append(tag)
            tracks.append(track)
            break
        
    
    for probe in range(evt_tree.recoMuon_size):
        ####################Probe Denominator Cuts#########################################
        probe_eta = evt_tree.recoMuon_etaSt2[probe]
        probe_phi = evt_tree.recoMuon_phiSt2[probe]
        isStation2 = True
        if probe_eta < -99 or probe_phi < -99: 
            probe_eta = evt_tree.recoMuon_etaSt1[probe]
            probe_phi = evt_tree.recoMuon_phiSt1[probe]
            isStation2 = False
        if probe_eta < -99 or probe_phi < -99: continue 

        probe_phi_deg = probe_phi * 180. / 3.14159
        probe_pt = evt_tree.recoMuon_pt[probe]
        probe_pt = min(probe_pt, max_pt)
        
        #require specific tightness
        if not REQ_TIGHT and (not evt_tree.recoMuon_isMediumMuon[probe] or abs(evt_tree.recoMuon_dxy[probe]) > .2 or abs(evt_tree.recoMuon_dz[probe]) > .2): continue
        if REQ_TIGHT and not evt_tree.recoMuon_isTightMuon[probe]: continue
        if MEDIUM_ONLY and evt_tree.recoMuon_isTightMuon[probe]: continue

        #########Tag/Probe Isolation #####################
        matched_tag = -1
        for tag in tags:

            tag_eta = evt_tree.recoMuon_etaSt2[tag]
            tag_phi = evt_tree.recoMuon_phiSt2[tag]
            if tag_eta < -99 or tag_phi < -99: 
                tag_eta = evt_tree.recoMuon_etaSt1[tag]
                tag_phi = evt_tree.recoMuon_phiSt1[tag]

            if probe == tag: continue
            if run >= 356798 and evt_tree.emtfUnpTrack_bx[track] != -1: continue
            elif run < 356798 and evt_tree.emtfUnpTrack_bx[track] != 0: continue

            #require probe has a distinct tag, ensure probe is not the same muon as the tag muon
            if h.CalcDR(tag_eta, tag_phi, probe_eta, probe_phi) < 4*MAX_dR: continue

            matched_tag = tag
            break
        ##################################################

        if matched_tag < 0: continue
        

        #require probe muon to be in EMTF range
        if abs(probe_eta) < eta_min or abs(probe_eta) > eta_max: continue

        #Fill denominator only with muons in the efficiency plateau (pt > 26 GeV)
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
        

        #Match probe to a track
        best_track = -1
        best_dr = -1
        for track in range(evt_tree.emtfUnpTrack_size):

            track_eta = evt_tree.emtfUnpTrack_eta[track]
            track_phi = evt_tree.emtfUnpTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159
            

            new_dr = h.CalcDR( track_eta, track_phi, probe_eta, probe_phi)

            
            #Ensure track does not belong to tag
            if tracks[tags.index(matched_tag)] == track: continue
            if run >= 356798 and evt_tree.emtfUnpTrack_bx[track] != -1: continue
            elif run < 356798 and evt_tree.emtfUnpTrack_bx[track] != 0: continue
            #Track is sufficiently close to probe
            if new_dr > 2 * MAX_dR: continue
            
            #If first track or closest track, assign as the best track
            if best_track == -1 or new_dr < best_dr: 
                best_dr = new_dr
                best_track = track


        ##############TRIG(Numerator) CUTS###################

        #Ensure track exists
        if best_track < 0: continue
        track = best_track
        
        track_eta = evt_tree.emtfUnpTrack_eta[track]
        track_phi = evt_tree.emtfUnpTrack_phi[track]
        track_phi *= 3.14159 / 180.0
        mode = evt_tree.emtfUnpTrack_mode[track]

        ##Fill denoinator plots (no PT cut) for different Mode Cuts, then fill trigger plots if pt  > TRG_PT
        for i in range(len(modes)):
            if probe_pt > PRB_PT and evt_tree.emtfUnpTrack_mode[track] in modes[i]: 
                if probe_eta < -2.2:
                    phi_mode[0][i].Fill(probe_phi)
                    if evt_tree.emtfUnpTrack_pt[track] > TRG_PT:
                        phi_mode_trig[0][i].Fill(probe_phi)
                elif probe_eta > 2.2:
                    phi_mode[1][i].Fill(probe_phi)
                    if evt_tree.emtfUnpTrack_pt[track] > TRG_PT:
                        phi_mode_trig[1][i].Fill(probe_phi)
        
        ##Fill trigger plots (no PT cut) for different PT Cuts
        for i in range(len(cuts)):
            if evt_tree.emtfUnpTrack_pt[track] > cuts[i] and probe_pt > PRB_PT and mode in [11, 13, 14, 15]:
                
                if probe_eta < -2.2:
                    phi_cut_trig[0][i].Fill(probe_phi)
                    
                elif probe_eta > 2.2:
                    phi_cut_trig[1][i].Fill(probe_phi)


        #Fill trigger plots for differently charge muons, standard mode and PT cuts
        if evt_tree.emtfUnpTrack_pt[track] > TRG_PT and probe_pt > PRB_PT and mode in [11, 13, 14, 15]:
            if probe_eta > 2.2:
                if evt_tree.recoMuon_charge[probe] < 0:
                    h_phi_neg_trig.Fill(probe_phi)
                elif evt_tree.recoMuon_charge[probe] > 0:
                    h_phi_pos_trig.Fill(probe_phi)


            #Station-Station DPhi plots (we do this with CSC hits or with phi to see which geometry is more accurate)
            #Find which station-station transsitions exist based on the mode
            for id in d_phi_mode_map[mode]:
                
                #Get the two station hits
                first_station = d_phi_index_map[id][0]
                second_station = d_phi_index_map[id][1]

                hitref_1 = eval("evt_tree.emtfUnpTrack_hitref%d[track]" % (first_station))
                hitref_2 = eval("evt_tree.emtfUnpTrack_hitref%d[track]" % (second_station))
                
                # #Match hits to segments
                # seg1 = evt_tree.emtfUnpHit_match_iSeg[hitref_1]
                # seg2 = evt_tree.emtfUnpHit_match_iSeg[hitref_2]


                # #Get global segment phi from rectangular coordinates
                # globy = evt_tree.cscSegment_globY[seg1]
                # globx = evt_tree.cscSegment_globX[seg1]
                # csc_phi1 = math.atan(globy / globx)
                # if globx < 0:
                #     csc_phi1 += 3.14159
                # if csc_phi1 > 3.14159: csc_phi1 -= 2 * 3.14159
                # elif csc_phi1 < -3.14159: csc_phi1 += 2 * 3.14159 

                # #Do same with second segment
                # globy = evt_tree.cscSegment_globY[seg2]
                # globx = evt_tree.cscSegment_globX[seg2]
                # csc_phi2 = math.atan(globy / globx)
                # if globx < 0:
                #     csc_phi2 += 3.14159
                # if csc_phi2 > 3.14159: csc_phi2 -= 2 * 3.14159
                # elif csc_phi2 < -3.14159: csc_phi2 += 2 * 3.14159 

                # #Find difference in phi's between CSC hits
                # d_phi = csc_phi2 - csc_phi1
                # d_phi *= 180./3.14159

                # if d_phi < -180: d_phi += 360.
                # elif d_phi > 180: d_phi -= 360.

                # d_phi *= 60 #60 parts per degree


                #Comment this line out to use delta phi for CSC segments
                d_phi = evt_tree.emtfUnpHit_emtf_phi[hitref_2] - evt_tree.emtfUnpHit_emtf_phi[hitref_1]
                charge = evt_tree.emtfUnpTrack_q[track] == 1
                #Fill d_phi by phi, indexed station-station id, sector, and endcap
                d_phi_plots[evt_tree.emtfUnpTrack_endcap[track] == 1][id][0][charge].Fill(d_phi)
                
                d_phi_plots[evt_tree.emtfUnpTrack_endcap[track] == 1][id][evt_tree.emtfUnpTrack_sector[track]][charge].Fill(d_phi)

                

        #Loop through 4 stations to get Dphi between hit and matched segment 
        #(NOTE THAT THIS IS ONLY FOR TAGGED EVENTS, SEE alignment_study.py FOR BETTER STATISTICS)
        # for n in range(4):

            # #Get track hit and matched segment
            # hit_ref = eval("evt_tree.emtfUnpTrack_hitref%d[track]" % (n + 1))
            # if hit_ref < 0: continue
            # iSeg = evt_tree.emtfUnpHit_match_iSeg[hit_ref]
            # if iSeg < 0: continue

            # #Get hit information
            # hit_theta = evt_tree.emtfUnpHit_emtf_theta[hit_ref]
            # hit_theta *= 3.14159 / 180.0
            # hit_eta = -math.log(math.tan(hit_theta/2.0))
            # hit_sector = evt_tree.emtfUnpHit_sector[hit_ref]

            # #Convert hit phi to global coordinates
            # hit_phi_loc = evt_tree.emtfUnpHit_emtf_phi[hit_ref] / 60.0 - 22.
            # hit_phi_glob = hit_phi_loc + 15. + (60. * (hit_sector - 1))
            # if hit_phi_glob > 180.0: hit_phi_glob -= 360.
            # hit_phi = hit_phi_glob * 3.14159/180.
            

            # #Get CSC coordinates
            # csc_eta = evt_tree.cscSegment_eta[iSeg]
            # csc_phi = evt_tree.cscSegment_phi[iSeg] * 3.14159

            # globy = evt_tree.cscSegment_globY[iSeg]
            # globx = evt_tree.cscSegment_globX[iSeg]
            # csc_phi = math.atan(globy / globx)

            # if globx < 0:
            #     csc_phi += 3.14159
            
            # if csc_phi > 3.14159: csc_phi -= 2 * 3.14159
            # elif csc_phi < -3.14159: csc_phi += 2 * 3.14159 


            # #Get phi differerence betweeen CSC and hit
            # dphi = csc_phi - hit_phi
            
            # if dphi > 3.14159: dphi -= 2 * 3.14159
            # elif dphi < -3.14159: dphi += 2 * 3.14159
            # dphi *= 180. / 3.14159
            # #print('Glob x ' + str(globx) + 'Glob y  ' + str(globy) + ' ' + str(dphi) + 'Eta ' + str(track_eta) + 'Phi ' + str(track_hi))
            # #print('Event' + str(event) + ' Track ' + str(track) + ' Hit ' + str(hit_ref) + ' ' + 'Probe ' + str(probe) + ' Eta ' + str(track_eta) + ' Phi ' + str(track_phi) + ' Mode ' + str(mode))
            # ring = (evt_tree.cscSegment_ring[iSeg] - 1) % 3
            
            # #Plot by sattion, and ring wrt chamber, only 18 chambers when station > 1 and ring == 2
            # if n > 0 and ring == 0:
            #     phi_diff_plot[evt_tree.cscSegment_endcap[iSeg] == 1][n][ring].Fill(2 * evt_tree.cscSegment_chamber[iSeg] - 1, dphi)
            #     phi_diff_plot[evt_tree.cscSegment_endcap[iSeg] == 1][n][ring].Fill(2 * evt_tree.cscSegment_chamber[iSeg], dphi)
            # else:
            #     phi_diff_plot[evt_tree.cscSegment_endcap[iSeg] == 1][n][ring].Fill(evt_tree.cscSegment_chamber[iSeg], dphi)
            
        

#Statistical error graph of pt efficiencies
err_phi_mode = [[], []]

#Write Efficiency vs. phi with indexed mode cuts
for i in range(len(modes)):
    for j in range(2):

        err_phi_mode[j].append(TGraphAsymmErrors(phi_mode_trig[j][i], phi_mode[j][i]))
        err_phi_mode[j][i].SetName('h_phi_EMTF_%s_mode_err_%d' % (str(modes[i]), j))
        err_phi_mode[j][i].GetXaxis().SetTitle("Phi");
        err_phi_mode[j][i].GetYaxis().SetTitle("L1T Efficiency");
        err_phi_mode[j][i].GetYaxis().SetRangeUser(0.00001,1.2);

        canvas = TCanvas(err_phi_mode[j][i].GetName() , err_phi_mode[j][i].GetName(), 700,700)
        err_phi_mode[j][i].Draw("A P C")
    
        err_phi_mode[j][i].SetLineWidth(1)
        err_phi_mode[j][i].SetLineColor(1)

        cms_label =cms_latex()
        header = head()

        gStyle.SetLegendBorderSize(0)
        gStyle.SetLegendTextSize(0.018);

        leg =TLegend(0.4,0.8,0.88,0.88);
        if j == 0:
            leg.AddEntry(err_phi_mode[j][i], "L1 pT (prompt) > 22 GeV, eta < -2.2")
        if j == 1:
            leg.AddEntry(err_phi_mode[j][i], "L1 pT (prompt) > 22 GeV, eta > 2.2")
        leg.Draw("same P")
        
        #Write PDF's for error efficiency plots
        canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/" + err_phi_mode[j][i].GetName() + "unp.pdf")
        err_phi_mode[j][i].Write()

        phi_mode_trig[j][i].Write()
        phi_mode[j][i].Write()
        phi_mode_trig[j][i].Divide(phi_mode[j][i])
        phi_mode_trig[j][i].SetName('phi_mode_eff_unp')
        phi_mode_trig[j][i].Write()



#Write Efficiency vs. Phi indexed by trigger PT cuts
err_phi_cut = [[], []]
for i in range(len(cuts)):
    for j in range(2):
        err_phi_cut[j].append(TGraphAsymmErrors(phi_cut_trig[j][i], phi_cut[j][i]))
        err_phi_cut[j][i].GetYaxis().SetRangeUser(0, 1)
        err_phi_cut[j][i].SetName('h_phi_EMTF_%s_cut_err_%d' % (str(cuts[i]), j))
        err_phi_cut[j][i].GetXaxis().SetTitle("Phi");
        err_phi_cut[j][i].GetYaxis().SetTitle("L1T Efficiency");
        err_phi_cut[j][i].GetYaxis().SetRangeUser(0.00001,1.2);

        canvas = TCanvas(err_phi_cut[j][i].GetName() , err_phi_cut[j][i].GetName(), 700,700)
        err_phi_cut[j][i].Draw("A P C")


    
        err_phi_cut[j][i].SetLineWidth(1)
        err_phi_cut[j][i].SetLineColor(1)

        cms_label =cms_latex()
        header = head()

        gStyle.SetLegendBorderSize(0)
        gStyle.SetLegendTextSize(0.018);

        leg =TLegend(0.4,0.8,0.88,0.88);
        if j == 0:
            leg.AddEntry(err_phi_cut[j][i], "L1 pT (prompt) > 22 GeV, eta < -2.2")
        if j == 1:
            leg.AddEntry(err_phi_cut[j][i], "L1 pT (prompt) > 22 GeV, eta > 2.2")
        leg.Draw("same P")
        

        #Save error plots as PDFs
        canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/" + err_phi_cut[j][i].GetName() + "unp.pdf")
        err_phi_cut[j][i].Write()



        phi_cut_trig[j][i].Write() 
        phi_cut[j][i].Write()
        phi_cut_trig[j][i].Divide(phi_cut[j][i])
        phi_cut_trig[j][i].SetName('phi_cut_eff_unp')
        phi_cut_trig[j][i].Write()    




#Write efficiency plots in high-eta for negative muons
err_phi_neg = TGraphAsymmErrors(h_phi_neg_trig, h_phi_neg)
h_phi_neg_trig.Write()
h_phi_neg.Write()
err_phi_neg.GetYaxis().SetRangeUser(0, 1)
err_phi_neg.SetName('h_phi_EMTF_neg_err')
err_phi_neg.GetXaxis().SetTitle("Phi");
err_phi_neg.GetYaxis().SetTitle("L1T Efficiency");
err_phi_neg.GetYaxis().SetRangeUser(0.00001,1.2);

canvas = TCanvas(err_phi_neg.GetName() , err_phi_neg.GetName(), 700,700)
err_phi_neg.Draw("A P C")

err_phi_neg.SetLineWidth(1)
err_phi_neg.SetLineColor(1)

cms_label =cms_latex()
header = head()

gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018);

leg =TLegend(0.4,0.8,0.88,0.88);
leg.AddEntry(err_phi_neg, "L1 pT (prompt) > 22 GeV, eta > 2.2, Negative Muon")
leg.Draw("same P")
canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/" + err_phi_neg.GetName() + "unp.pdf")




#Write efficiency plots in high-eta for negative muons
err_phi_pos = TGraphAsymmErrors(h_phi_pos_trig, h_phi_pos)
h_phi_pos_trig.Write()
h_phi_pos.Write()
err_phi_pos.SetName('h_phi_EMTF_pos_err')
err_phi_pos.GetXaxis().SetTitle("Phi");
err_phi_pos.GetYaxis().SetTitle("L1T Efficiency");
err_phi_pos.GetYaxis().SetRangeUser(0.00001,1.2);

canvas = TCanvas(err_phi_pos.GetName() , err_phi_pos.GetName(), 700,700)
err_phi_pos.Draw("A P C")

err_phi_pos.SetLineWidth(1)
err_phi_pos.SetLineColor(1)

cms_label =cms_latex()
header = head()

gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018);

leg =TLegend(0.4,0.8,0.88,0.88);
leg.AddEntry(err_phi_pos, "L1 pT (prompt) > 22 GeV, eta > 2.2, Positive Muon")
leg.Draw("same P")
canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/" + err_phi_pos.GetName() + "unp.pdf")



#Write station-station phi differences wrt phi by sector
for endcap in range(2):
    for d in range(6):
        for sector in range(7):
            for charge in range(2):
                d_phi_plots[endcap][d][sector][charge].Write()


#Write CSC segment and hit phi differences wrt phi
for endcap in range(2):
    for station in range(4):
        for ring in range(3):
            phi_diff_plot[endcap][station][ring].Write()
    

del out_file
