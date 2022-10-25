from cmath import tan
from re import L
import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
import Helper as h
import argparse

''' This script exists to access all hits and tracks in events for the purpose of analyzing potential Run 3 misalignments'''

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


parser = argparse.ArgumentParser()
parser.add_argument("-n", "--num_jobs", required=False)
parser.add_argument("-i", "--index", required = False)
args = parser.parse_args()



## Configuration settings
MAX_FILE =  -1   ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event

#Specify which datasets to use
IDEAL = False
ZERO_BIAS = False
RUN3 = True



MAX_dR = .1 #Proximity for track matching
s6_chambers = [33, 34, 35, 36, 1, 2] #projections from chambers to sector 6 require some fudging since they aren't contiguous
evt_tree  = TChain('EMTFNtuple/tree') #event tree for analyzing event data

#parameters by which plots are indexed
pt_ranges = [(0, 5), (5, 10), (10, 20), (20, 50), (0, 50), (0, 20)]
eta_ranges = [(1.2, 1.6), (1.6, 2.1,), (2.1, 2.4)]
modes = [range(16), [11, 13, 14, 15], [7, 9, 10, 11, 13, 14, 15], [7], [9], [10], [11], [13], [14], [15]]


#arrays of plots
reso_plots = [] #Plots pt-resolution with respect to eta and phi
reso_plots_1d = [] #plots pt-resolution purely with respect to phi in the eta > 2.2 range
sum_plots = [] #sums the number of hits in a given eta-phi bin for normalization
sum_plots_1d = [] #sums the number of hits in a given phi bin for normalization

#histogram of pt-resolutions, indexed by pt-range, eta-region, endcap, and modes of tracks
reso_histo = np.zeros((len(pt_ranges) + 1, len(eta_ranges), 2, len(modes)), dtype=object)

#reso plots are indexed by pt_ranges
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

    #resolution gaussian hisograms are also indexed by eta ranges, endcaps, and track modes
    for j, bound_eta in enumerate(eta_ranges):
        for endcap in range(2):
            for k in range(len(modes)):
                reso_histo[i][j][endcap][k] = TH1D('pt_resolution_histo_e%d_%d_%d_%d_%d_%s_%s' % 
                    (endcap,
                    bound[0],
                    bound[1],
                    bound_eta[0],
                    bound_eta[1],
                    str(modes[k]),
                    "ideal" if IDEAL else ""
                    ), '', 60, -3.01, 3.01)

                #Last element of the array includes all PT's, only allow this to be filled for one pt-range iteration
                if i == 0:
                    reso_histo[len(pt_ranges)][j][endcap][k] = TH1D('pt_resolution_histo_e%d_anypt_%d_%d_%s_%s' % 
                        (endcap,
                        bound_eta[0],
                        bound_eta[1],
                        str(modes[k]),
                        "ideal" if IDEAL else ""
                        ), '', 60, -3.01, 3.01)


#encap, station, ring
phi_diff_plot = np.zeros(shape=(2, 4, 3), dtype=object) #plot csc phi - hit phi for different chambers (indexed per endcap, station, ring)
for endcap in range(2):
    for station in range(4):
        for ring in range(3):
                phi_diff_plot[endcap][station][ring] = TH2D('phi_diff_EMTF_%d_%d_%d' % (endcap, station + 1, ring + 1), '', 36, 1, 37, 100, -.5, .5)
                phi_diff_plot[endcap][station][ring].GetYaxis().SetTitle("LCT #phi - Hit #phi")
                phi_diff_plot[endcap][station][ring].GetXaxis().SetTitle("Chamber")
                phi_diff_plot[endcap][station][ring].SetTitle("#Delta#phi Endcap %d, Station %d, Ring %d" %
                (1 if endcap else -1, 
                station + 1,
                ring + 1))






if args.num_jobs:
    INDEX = int(args.index)
    NUM_JOBS = int(args.num_jobs)

#the gaussian mean of phi_diff (per chamber/sector), and find difference in these chamber/sector means between stations
#this will give the total misaligned delta phi for station-station transitions
station_station_plot_chamber = np.zeros(shape=(2, 4, 4, 3), dtype=object) #gaussian mean projected per chamber
station_station_plot = np.zeros(shape=(2, 4, 4, 3), dtype=object) #gaussian mean projected per sector

track_dphi = np.zeros(shape=(2, 4, 4, 3, 2), dtype=object) #gaussian mean projected per sector
for endcap in range(2):
    for station1 in range(4):
        for station2 in range(station1 + 1, 4):
            for ring in range(3):
                station_station_plot[endcap][station1][station2][ring][charge] = TH1D('sts_e%d_(%d,%d)_r%d' % (endcap, station1 + 1, station2 + 1, ring + 1), '', 6, 1, 7)
                station_station_plot[endcap][station1][station2][ring].GetYaxis().SetTitle('#Delta#phi_{%d} - #Delta#phi_{%d}' % (station2 + 1, station1 + 1))
                station_station_plot[endcap][station1][station2][ring].GetXaxis().SetTitle("Sector")
                station_station_plot[endcap][station1][station2][ring].SetTitle("#Delta#phi_{%d%d}, Endcap %d, Ring %d" %
                (station1 + 1, 
                station2 + 1,
                1 if endcap else -1, 
                ring + 1))

                station_station_plot_chamber[endcap][station1][station2][ring] = TH1D('sts_chamber_e%d_(%d,%d)_r%d' % (endcap, station1 + 1, station2 + 1, ring + 1), '', 36, 1, 37)
                station_station_plot_chamber[endcap][station1][station2][ring].GetYaxis().SetTitle('#Delta#phi_{%d} - #Delta#phi_{%d}' % (station2 + 1, station1 + 1))
                station_station_plot_chamber[endcap][station1][station2][ring].GetXaxis().SetTitle("Chamber")
                station_station_plot_chamber[endcap][station1][station2][ring].SetTitle("#Delta#phi_{%d%d}, Endcap %d, Ring %d" %
                (station1 + 1, 
                station2 + 1,
                1 if endcap else -1, 
                ring + 1))
            


if args.num_jobs:
  try: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/alignment_study_unpacked%d.root" % (INDEX), 'create')
  except: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/alignment_study_unpacked%d.root" % (INDEX), 'recreate')
else:
  out_file = TFile('plots/alignment/alignment_study_unpacked_single.root', 'recreate')

#Determine whether to use ZeroBias dataset and ideal geometry look-up tables

if ZERO_BIAS:
    if not IDEAL:
        folder = "/eos/user/n/nhurley/ZeroBias/EMTFNtuple_Run3_ZeroBias_data_13p6TeV_Run2Alignment/220918_123516/0000/"

    if IDEAL: 
        folder = "/eos/user/n/nhurley/ZeroBias/EMTFNtuple_Run3_ZeroBias_data_13p6TeV_idealAlignment/220918_123144/0000/"
elif not ZERO_BIAS:
    
    if IDEAL: 
        folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_SingleMuon_data_13p6TeV_idealAlignment/220902_142649/0000/"
    if RUN3:
        folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run3Alignment_2022C_v5/220925_185603/0000/"
  

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
    if nFiles == MAX_FILE: break

#Go event by event
for event in range(evt_tree.GetEntries()):
    evt_tree.GetEntry(event)

    if event % PRT_EVT == 0:
        if INDEX:
            print('alignment_study_unpacked.py: Processing Job #%d, Event #%d' % (INDEX, event))
        else: print('alignment_study_unpacked.py: Processing Event #%d' % (event))

    #Analyze all hits for phi_diff plots
    # for hit_ref in range(evt_tree.emtfUnpHit_size):
        
    #     #Exclude neighbor hits
    #     if evt_tree.emtfUnpHit_neighbor[hit_ref]: continue

    #     #Match to a CSC Segment
    #     iSeg = evt_tree.emtfUnpHit_match_iSeg[hit_ref]
    #     if iSeg < 0: continue

    #     #Calculate hit coordinates
    #     hit_theta = evt_tree.emtfUnpHit_emtf_theta[hit_ref]
    #     hit_theta *= 3.14159 / 180.0
    #     hit_eta = -math.log(math.tan(hit_theta/2.0))
        
    #     #Only analyze hits in high-eta region
    #     #if hit_eta < 2.2: continue

    #     #Convert hit-phi to global phi in radians
    #     hit_sector = evt_tree.emtfUnpHit_sector[hit_ref]
    #     hit_phi_loc = evt_tree.emtfUnpHit_emtf_phi[hit_ref] / 60.0 - 22.
    #     hit_phi_glob = hit_phi_loc + 15. + (60. * (hit_sector - 1))
    #     if hit_phi_glob > 180.0: hit_phi_glob -= 360.
    #     hit_phi = hit_phi_glob * 3.14159/180.
        

    #     #Convert CSC rectangular coordinates to phi
    #     csc_eta = evt_tree.cscSegment_eta[iSeg]
    #     csc_phi = evt_tree.cscSegment_phi[iSeg] * 3.14159
    #     globy = evt_tree.cscSegment_globY[iSeg]
    #     globx = evt_tree.cscSegment_globX[iSeg]
    #     csc_phi = math.atan(globy / globx)
    #     if globx < 0:
    #         csc_phi += 3.14159
    #     if csc_phi > 3.14159: csc_phi -= 2 * 3.14159
    #     elif csc_phi < -3.14159: csc_phi += 2 * 3.14159 

    #     #Get differentce between CSC (accurate) and hit (fast) phi's
    #     dphi = csc_phi - hit_phi
    #     if dphi > 3.14159: dphi -= 2 * 3.14159
    #     elif dphi < -3.14159: dphi += 2 * 3.14159
    #     dphi *= 180. / 3.14159

    #     #Convert ring / station to indicies
    #     ring = (evt_tree.cscSegment_ring[iSeg] - 1) % 3
    #     station = evt_tree.cscSegment_station[iSeg] - 1

    #     #Station 2/3/4 Ring 1 have 18 chambers instead of 36
    #     if station > 0 and ring == 0:
    #         phi_diff_plot[evt_tree.cscSegment_endcap[iSeg] == 1][station][ring].Fill(2 * evt_tree.cscSegment_chamber[iSeg] - 1, dphi)
    #         phi_diff_plot[evt_tree.cscSegment_endcap[iSeg] == 1][station][ring].Fill(2 * evt_tree.cscSegment_chamber[iSeg], dphi)
    #     else:
    #         phi_diff_plot[evt_tree.cscSegment_endcap[iSeg] == 1][station][ring].Fill(evt_tree.cscSegment_chamber[iSeg], dphi)

    #Match all reconstructed muons to tracks (without tagging and cutting them), called probe as a holdover from tag&probe studies
    for probe in range(evt_tree.recoMuon_size):
        #Cut by reconstructed muon quality
        if not evt_tree.recoMuon_isMediumMuon[probe]: continue
        probe_eta = evt_tree.recoMuon_etaSt2[probe]
        probe_phi = evt_tree.recoMuon_phiSt2[probe]

        #make sure reco muon eta is in station 2, if not station 1
        if probe_eta < -99 or probe_phi < -99: 
            probe_eta = evt_tree.recoMuon_etaSt1[probe]
            probe_phi = evt_tree.recoMuon_phiSt1[probe]
        if probe_eta < -99 or probe_phi < -99: continue 
        probe_phi_deg = probe_phi * 180. / 3.14159
        probe_pt = evt_tree.recoMuon_pt[probe]

        #Find track with the lowest distance to the reco-muon
        best_track = -1
        best_dr = -1
        for track in range(evt_tree.emtfUnpTrack_size):

            track_eta = evt_tree.emtfUnpTrack_eta[track]
            track_phi = evt_tree.emtfUnpTrack_phi[track]
            track_phi *= 3.14159 / 180.0
            if (track_phi > 3.14159): track_phi -= 2*3.14159
            elif (track_phi < -3.14159): track_phi += 2*3.14159

            
            ##########This section can be used to check alignment##################################
            # for unp_track in range(evt_tree.emtfUnpTrack_size):
            #     unp_track_eta = evt_tree.emtfUnpTrack_eta[unp_track]
            #     unp_track_phi = evt_tree.emtfUnpTrack_phi[unp_track]
            #     unp_track_phi *= 3.14159 / 180.0
            #     if (unp_track_phi > 3.14159): unp_track_phi -= 2*3.14159
            #     elif (unp_track_phi < -3.14159): unp_track_phi += 2*3.14159

            #     if h.CalcDR( track_eta, track_phi, unp_track_eta, unp_track_phi) < .2:
            #         if unp_track_phi != track_phi:
            #             print("Unpacked phi: " + str(unp_track_phi) + ", Re-emulated phi: " + str(track_phi))
            #         break

            new_dr = h.CalcDR( track_eta, track_phi, probe_eta, probe_phi)

            #If best track is unassigned or new track is closer to muon than the best track, then re-assign best-track
            if best_track == -1 or new_dr < best_dr:
                best_dr = new_dr
                best_track = track

        #if best track exceeds maximum distance from reco-muon, or isn't found, then skip this reco-muon
        if best_dr > 2 * MAX_dR or best_track < 0: continue
        track = best_track

        #Estimate of how well trigger-PT predicts reconstructed PT
        #Values < -1 are due to charge misidentification
        precision = (evt_tree.emtfUnpTrack_q[track] * evt_tree.recoMuon_charge[probe] * evt_tree.emtfUnpTrack_pt[track] - probe_pt) / probe_pt

        #Index PT-resolution plots by pt_ranges
        for i, bound in enumerate(pt_ranges):
            if bound[0] < probe_pt < bound[1]: 
                reso_plots[i].Fill(probe_eta, probe_phi, precision * precision)
                sum_plots[i].Fill(probe_eta, probe_phi)

                #1D reso plots only include high-eta region
                if probe_eta > 2.2:
                    reso_plots_1d[i].Fill(probe_phi, precision * precision)
                    sum_plots_1d[i].Fill(probe_phi)


                #PT gaussian histrograms are also indexed by mode, eta region, and endcap
                for j, bound_eta in enumerate(eta_ranges):
                    if bound_eta[0] < abs(probe_eta) < bound_eta[1]:
                        for k, mode in enumerate(modes):
                            if evt_tree.emtfUnpTrack_mode[track] in mode:
                                reso_histo[i][j][probe_eta > 0][k].Fill(precision)
                                
        for j, bound_eta in enumerate(eta_ranges):
                if bound_eta[0] < abs(probe_eta) < bound_eta[1]:
                    for k, mode in enumerate(modes):
                        if evt_tree.emtfUnpTrack_mode[track] in mode:
                            reso_histo[len(pt_ranges)][j][probe_eta > 0][k].Fill(precision)
                

#Write plots to file
for endcap in range(2):
    for station in range(4):
        for ring in range(3):
            phi_diff_plot[endcap][station][ring].Write()




gStyle.SetOptFit(1)
gauss = np.zeros(shape=(2, 4, 3, 6), dtype=object) #Array Gaussian means by sector
gauss_chamber = np.zeros(shape=(2, 4, 3, 36)) #Array Gaussian means by chamber
for endcap in range(2):
    for station in range(4):
        for ring in range(3):

            #Get 2D DPhi plot
            plot = phi_diff_plot[endcap][station][ring]

            #Project each chamber to a 1D gaussian and find its mean to get the mean DPhi
            for chamber in range(36):
                gauss_fit = plot.ProjectionY(
                            'proj_%d_%d_%d_%d' % (endcap, station, ring, chamber),
                            chamber + 1,
                            chamber + 1
                            ).Fit(formula='gaus', option="Q S")
                #If fit exists, assign mean to array
                if int(gauss_fit) != -1:
                    gauss_chamber[endcap][station][ring][chamber] = gauss_fit.Parameters()[1]


                
            #Project each sector to a 1D gaussian and save its mean
            for sector in range(6):
                
                #Sector 6 has non-contiguous plots, so we improvise
                if sector == 5:
                    #Plot to be fitted to
                    sector_six = TH1D('sector_six_%d_%d_%d' % (endcap, station, ring), '', 100, -.5, .5)

                    #Project each chamber to a 1D plot in sector 6, and add it to the total sector plot
                    for chamber in s6_chambers:
                        sector_six.Add(plot.ProjectionY(
                            'proj_chamber_%d_%d_%d_%d_%d' % (endcap, station, ring, sector, chamber),
                            chamber, chamber))
                    #Find gaussian fit of the sum of the plots in each sector-six chamber
                    gauss_fit = sector_six.Fit(formula='gaus', option='Q S')


                #All other sectors can have all 6 chambers projected with a single line
                else:
                    gauss_fit = plot.ProjectionY(
                            'proj_%d_%d_%d_%d' % (endcap, station, ring, sector),
                            3 + 6*sector,
                            (2 + 6 * (sector + 1)) % 36
                            ).Fit(formula='gaus', option="Q S")
                #Save gaussian fit if it exists
                if int(gauss_fit) == -1:
                    gauss[endcap][station][ring][sector] = 0
                else:
                    gauss[endcap][station][ring][sector] = gauss_fit.Parameters()[1]


#Create station to station plots with these gaussian means
for endcap in range(2):
    for station1 in range(4):
        for station2 in range(station1 + 1, 4):
            for ring in range(3):
                
                #Find chamber means at different stations
                for chamber in range(36):
                    mean1 = gauss_chamber[endcap][station1][ring][chamber]
                    mean2 = gauss_chamber[endcap][station2][ring][chamber]
                    
                    #Default to 0 if a chamber has no hits
                    if mean1 == 0 or mean2 == 0:
                        difference = 0
                    #Take difference between subsequent station and previous station
                    else: difference = gauss_chamber[endcap][station2][ring][chamber] - gauss_chamber[endcap][station1][ring][chamber]

                    #Fill plot by chamber in x-axis
                    station_station_plot_chamber[endcap][station1][station2][ring].Fill(chamber + 1, difference)
                
                #Fill plot by sector in x-axis
                for sector in range(6):
                    difference = gauss[endcap][station2][ring][sector] - gauss[endcap][station1][ring][sector]
                    station_station_plot[endcap][station1][station2][ring].Fill(sector + 1, difference)


#Write the station-station plots to file
for endcap in range(2):
    for station1 in range(4):
        for station2 in range(station1 + 1, 4):
            for ring in range(3):
                    station_station_plot[endcap][station1][station2][ring].Write()
                    station_station_plot_chamber[endcap][station1][station2][ring].Write()



#Write 2D resolution plots to file and normalize
for i, plot in enumerate(reso_plots):
    plot.Write()
    sum_plots[i].Write()
    plot.Divide(sum_plots[i])
    plot.SetName(plot.GetName() + "_normal")
    plot.Write()
    plot.Draw("A P C")

#Write 1D resolution plots to file and normalize
for i, plot in enumerate(reso_plots_1d):
    plot.Divide(sum_plots_1d[i])
    plot.SetName(plot.GetName() + "_normal")
    plot.Write()
    plot.Draw("A P")
    
    #Write PT gaussian histogram plots to file
    for j in range(len(eta_ranges)):
        for endcap in range(2):
            for k in range(len(modes)):
                reso_histo[i][j][endcap][k].Write()
                reso_histo[i][j][endcap][k].Draw("A P C")
                #Last element of the array includes all PT's, only allow this to be filled for one pt-range iteration
                if i == 0:
                    reso_histo[len(pt_ranges)][j][endcap][k].Write()

del out_file