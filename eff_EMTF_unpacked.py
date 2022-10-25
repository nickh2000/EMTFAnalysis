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
RUN3 = True
RUN2 = False

REQ_BX0    = False  ## Require L1T muon to be in BX 0  ## Require a final uGMT candidate, not just a TF muon
REQ_TIGHT = False 
MEDIUM_ONLY = False ## Require that muons be of medium quality but not of tight quality
REQ_HLT = False

MAX_dR  = 0.1     ## Maximum dR for L1T-offline matching
TAG_ISO = 0.15   ## Maximum relative isolation for tag muon
TAG_PT  = 26    ## Minimum offline pT for tag muon
PRB_PT  = 26    ## Minimum offline pT for probe muon
TRG_PT  = 22    ## Minimum L1T pT for probe muon

eta_min = 1.24
eta_max = 2.40
modes = [7, 9, 10, 11, 13, 14, 15]

scale_pt_temp = [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 45, 60, 75, 100, 140, 160, 180, 200, 250, 300, 500, 1000]
scale_pt_temp_2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 55, 60]
scale_pt_2  = array('f', scale_pt_temp_2)
scale_pt  = array('f', scale_pt_temp)
max_pt = scale_pt_temp[len(scale_pt_temp) - 1] - 0.01

pt_ranges = [(0, 5), (5, 10), (10, 20), (20, 50), (0, 50)]
eta_ranges = [(1.2, 1.6), (1.6, 2.1,), (2.1, 2.4)]
reso_plots = []
reso_plots_1d = []
sum_plots = []
sum_plots_1d = []

reso_histo = np.zeros((len(pt_ranges), len(eta_ranges)), dtype=object)

for i, bound in enumerate(pt_ranges):
    reso_plots.append(TH2D('pt_resolution_2d_%d_%d%s' % 
    (bound[0],
    bound[1],
    "_ideal" if IDEAL else ""
    ), '', 50, -2.4, 2.4, 30, -3.14159, 3.14159))

    reso_plots_1d.append(TH1D('pt_resolution_1d_%d_%d%s' % 
    (bound[0],
    bound[1],
    "_ideal" if IDEAL else ""
    ), '', 4, -3.14159, 3.14159))

    sum_plots.append(TH2D('pt_hit_sum_2d_%d_%d%s' % 
    (bound[0],
    bound[1],
    "_ideal" if IDEAL else ""
    ), '', 50, -2.4, 2.4, 30, -3.14159, 3.14159))

    sum_plots_1d.append(TH1D('pt_hit_sum_1d_%d_%d%s' % 
    (bound[0],
    bound[1],
    "_ideal" if IDEAL else ""
    ), '', 4, -3.14159, 3.14159))

    for j, bound_eta in enumerate(eta_ranges):
        reso_histo[i][j] = TH1D('pt_resolution_histo_%d_%d_%d_%d%s' % 
            (bound[0],
            bound[1],
            bound_eta[0],
            bound_eta[1],
            "_ideal" if IDEAL else ""
            ), '', 60, -3.01, 3.01)

## HISTOGRAMS 
evt_tree  = TChain('EMTFNtuple/tree')
h_pt  = TH1D('h_pt_den_EMTF',  '', len(scale_pt_temp_2) - 1,  scale_pt_2)
h_pt_trg = TH1D('h_pt_num_EMTF',  '', len(scale_pt_temp_2)-1,  scale_pt_2)
h_pt_large = TH1D('h_pt_den_EMTF_large',  '', len(scale_pt_temp) - 1,  scale_pt)
h_pt_trg_large = TH1D('h_pt_num_EMTF_large',  '', len(scale_pt_temp) - 1,  scale_pt)


h_eta  = TH1D('h_eta_EMTF',  '', 50,  -2.5, 2.5)
h_eta_trg = TH1D('h_eta_trg_EMTF',  '', 50, -2.5, 2.5)



h_eta_pt = TH2D('h_eta_pt_EMTF', '', 50, -2.4, 2.4, 50, 0, 100)
h_eta_phi = TH2D('h_eta_phi_EMTF', '', 50, -2.4, 2.4, 30, -3.14159, 3.14159)
h_eta_isTight =TH2D('h_eta_tight_EMTF', '', 50, -2.4, 2.4, 2, 0, 2)
h_eta_isTight_denom =TH2D('h_eta_tight_denom_EMTF', '', 50, -2.4, 2.4, 2, 0, 2) #plot percent of tight / medium muons by eta, need to divide by totoal muons of each tightness
pos_end_mode = TH2D('pos_end_mode', '', 30, -3.14159, 3.14159, 16, 0, 16) ##plot region in eta > 2.3 with low efficiency

h_eta_qual = TH2D('h_eta_qual_EMTF', '', 50, -2.4, 2.4, 16, 0, 16) #plot the mode vs eta

h_eta_pt_trg = TH2D('h_eta_pt_trg_EMTF', '', 50, -2.4, 2.4, 50, 0, 100)
h_eta_phi_trg = TH2D('h_eta_phi_trg_EMTF', '', 50, -2.4, 2.4, 30, -3.14159, 3.14159)
h_eta_isTight_trg =TH2D('h_eta_tight_trg_EMTF', '', 50, -2.4, 2.4, 2, 0, 2)




#Total muons of each tightness gives a fraction of tight/medium muons in each eta bin
tight_0_tot =0 
tight_1_tot = 0

if args.num_jobs:
    INDEX = int(args.index)
    NUM_JOBS = int(args.num_jobs)

if args.num_jobs:
  try: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/eff_EMTF_unpacked%d.root" % (INDEX), 'recreate')
  except: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/eff_EMTF_unpacked%d.root" % (INDEX), 'recreate')
else:
  out_file = TFile('plots/efficiency/EMFTF_eff_unpacked_single.root', 'recreate')

if IDEAL: 
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_SingleMuon_data_13p6TeV_idealAlignment/220902_142649/0000/"
elif RUN3:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run3Alignment_2022C_v5/220925_185603/0000/"
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

    if event % PRT_EVT == 0:
        if args.index:
            print('eff_EMTF_unpacked.py: Processing Job #%d, Event #%d' % (INDEX, event))
        else: print('eff_EMTF_unpacked.py: Processing Event #%d' % (event))

    tags = []
    tracks = []

    run = evt_tree.eventInfo_run[0]

    if RUN2:
        if int(run) > 359356: continue
    elif RUN3:
        if int(run) <= 359: continue

    for tag in range(evt_tree.recoMuon_size):
        

        #Require muon tightness
        if not REQ_TIGHT and not evt_tree.recoMuon_isMediumMuon[tag]: continue
        if REQ_TIGHT and not evt_tree.recoMuon_isTightMuon[tag]: continue


        #Tag needs to be a certain quality to be a valid tag candidate
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

        for track in range(evt_tree.emtfUnpTrack_size):
            
            
            track_eta = evt_tree.emtfUnpTrack_eta[track]
            track_phi = evt_tree.emtfUnpTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159
            
            #Ensure L1 trigger is of sufficient quality and close to the reconstructed tag
            if (evt_tree.emtfUnpTrack_mode[track] < 11 or evt_tree.emtfUnpTrack_mode[track] == 12): continue
            if h.CalcDR( track_eta, track_phi, tag_eta, tag_phi ) > MAX_dR: continue
            if (evt_tree.emtfUnpTrack_pt[track] < TRG_PT): continue

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

            #require probe has a distinct tag
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


        #fill tightness and phi plots only with muons in the efficiency plateau, higher quality (pt > 26 GeV)
        if (probe_pt > PRB_PT):

            h_eta.Fill(probe_eta)
            h_eta_isTight.Fill(probe_eta, int(evt_tree.recoMuon_isTightMuon[probe]))
            h_eta_phi.Fill(probe_eta, probe_phi)
            
        

    #################################################################################
        
        best_track = -1
        best_dr = -1
        for track in range(evt_tree.emtfUnpTrack_size):

            track_eta = evt_tree.emtfUnpTrack_eta[track]
            track_phi = evt_tree.emtfUnpTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159

            new_dr = h.CalcDR( track_eta, track_phi, probe_eta, probe_phi)

            ##############TRIG(Numerator) CUTS###################
            if run >= 356798 and evt_tree.emtfUnpTrack_bx[track] != -1: continue
            elif run < 356798 and evt_tree.emtfUnpTrack_bx[track] != 0: continue
            if tracks[tags.index(matched_tag)] == track: continue

            if best_track == -1 or new_dr < best_dr:
                best_dr = new_dr
                best_track = track

        if best_dr > 2 * MAX_dR or best_track < 0: continue
        track = best_track


        precision = (evt_tree.emtfUnpTrack_q[track] * evt_tree.emtfUnpTrack_pt[track] - evt_tree.recoMuon_charge[probe] *probe_pt) / probe_pt

        for i, bound in enumerate(pt_ranges):
            if evt_tree.emtfUnpTrack_mode[track] != 15: continue
            if bound[0] < probe_pt < bound[1]: 
                reso_plots[i].Fill(probe_eta, probe_phi, precision * precision)
                
                sum_plots[i].Fill(probe_eta, probe_phi)

                if probe_eta > 2.2:
                    reso_plots_1d[i].Fill(probe_phi, precision * precision)
                    sum_plots_1d[i].Fill(probe_phi)

                for j, bound_eta in enumerate(eta_ranges):
                    if bound_eta[0] < probe_eta < bound_eta[1]:
                        reso_histo[i][j].Fill(precision)

        if probe_eta >= 2.3 and evt_tree.emtfUnpTrack_pt[track] > TRG_PT: pos_end_mode.Fill(probe_phi, evt_tree.emtfUnpTrack_mode[track]) #find incidence of trigger muons in erroneous high eta region

        if evt_tree.emtfUnpTrack_pt[track] > 22 and probe_pt > PRB_PT and evt_tree.emtfUnpTrack_mode[track] in modes: h_eta_phi_trg.Fill(probe_eta, probe_phi) #looser pt constraints (7GeV) and loser mode constraints to get more detail on phi plot
 
        #Tag and probe momentum/quality cuts
        if (evt_tree.emtfUnpTrack_mode[track] < 11 or evt_tree.emtfUnpTrack_mode[track] == 12): continue #standard single-muon quality cuts
        if evt_tree.emtfUnpTrack_pt[track] < TRG_PT - .01: continue #single L1 trigger pt cut
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
            
            d_phi_12 = list(evt_tree.emtfUnpTrack_ptLUT_deltaPh)[track][0]
            

        if REQ_BX0 and evt_tree.emtfUnpTrack_bx[track] != -1: continue #ensure muon is in BX0

        h_eta_qual.Fill(probe_eta, evt_tree.emtfUnpTrack_mode[track]) #fill mode plot without regard to mode cuts

print ('%s: %.1f +/- %.1f%%' % ( "EMTF", 100 * h_pt_trg.Integral() / h_pt.Integral(), 
                                    (100 * h_pt_trg.Integral() / h_pt.Integral()) * 
                                    math.sqrt(h_pt.Integral()) / h_pt.Integral() ))

#Statistical error graph of pt efficiencies

h_err_pt = TGraphAsymmErrors(h_pt_trg, h_pt)
h_err_pt.SetName('h_pt_trg_EMTF_err',)
h_err_pt.Write()

h_pt_trg.Write() 
h_pt_trg_large.Write()
h_pt_trg.Divide(h_pt)
h_pt_trg.SetName(h_pt_trg.GetName()+'_eff')

h_pt.Write()
h_pt_large.Write()


#2D plot is similar to above, includes eta axis
h_eta_pt_trg.Write() 
h_eta_pt_trg.Divide(h_eta_pt)
h_eta_pt_trg.SetName(h_eta_pt_trg.GetName()+'_eff')
h_eta_pt_trg.Write()



#Plot efficiency of eta with respect to phi (only in 26GeV region)
h_eta_phi_trg.Write() 
h_eta_phi_trg.Divide(h_eta_phi)
h_eta_phi_trg.SetName(h_eta_phi_trg.GetName()+'_eff')

canvas = TCanvas(h_eta_phi_trg.GetName() , h_eta_phi_trg.GetName(), 700,700)
h_eta_phi_trg.Draw("colz")

h_eta_phi_trg.SetLineWidth(1)
h_eta_phi_trg.SetLineColor(1)

cms_label =cms_latex()
header = head()

gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018)
gStyle.SetOptStat(0)

if IDEAL:
    canvas.SaveAs("plots/ideal_" + h_eta_phi_trg.GetName() + ".pdf")
elif RUN3:
    canvas.SaveAs("plots/run3_" + h_eta_phi_trg.GetName() + ".pdf")
else:
    canvas.SaveAs("plots/" + h_eta_phi_trg.GetName() + ".pdf")

h_eta_phi_trg.Write()

h_err_eta = TGraphAsymmErrors(h_eta_trg, h_eta)
h_err_eta.SetName('h_eta_trg_EMTF_err',)
h_err_eta.SetName('h_eta_EMTF_err')
h_err_eta.GetXaxis().SetTitle("Eta")
h_err_eta.GetYaxis().SetTitle("L1T Efficiency")
h_err_eta.GetYaxis().SetRangeUser(0.00001,1.2)
h_err_eta.GetXaxis().SetRangeUser(-2.4,2.4)

canvas = TCanvas(h_err_eta.GetName() , h_err_eta.GetName(), 700,700)
h_err_eta.Draw("A P C")

h_err_eta.SetLineWidth(1)
h_err_eta.SetLineColor(1)

cms_label =cms_latex()
header = head()

gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018)

leg =TLegend(0.4,0.8,0.88,0.88)

leg.AddEntry(h_err_eta, "L1 pT (prompt) > 22 GeV")

canvas.SaveAs("plots/" + h_err_eta.GetName() + ".pdf")
#err_phi_mode.Write()

h_eta_trg.Write() 
h_eta_trg.Divide(h_eta)
h_eta_trg.SetName(h_eta_trg.GetName()+'_eff')
h_eta_trg.Write()

#Plot number/efficiency of tight/medium muuons
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

#plot mode vs eta
h_eta_qual.Write()

#plot mode vs phi in the high positive eta region
pos_end_mode.Write()

for i, plot in enumerate(reso_plots):
    plot.Write()
    sum_plots[i].Write()
    plot.Divide(sum_plots[i])
    plot.SetName(plot.GetName() + "_normal")
    plot.Write()
    plot.Draw("A P C")

for i, plot in enumerate(reso_plots_1d):
    plot.Divide(sum_plots_1d[i])
    plot.SetName(plot.GetName() + "_normal")
    plot.Write()
    plot.Draw("A P")

    for j in range(len(eta_ranges)):
        reso_histo[i][j].Write()
        plot.Draw("A P C")

del out_file