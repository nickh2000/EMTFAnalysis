from cmath import tan
from re import L
import sys
import math
from subprocess import Popen,PIPE
from tkinter.tix import MAX
from ROOT import *
import numpy as np
from array import *
import Helper as h
import os
import argparse


''' This script exists to analyze the rates of EMTF tracks
Specifically focuses on findings the RPC hits persist for two BX's in certain chambers
'''
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--num_jobs", required=False)
parser.add_argument("-i", "--index", required = False)
args = parser.parse_args()

if args.index:
  INDEX = int(args.index)
  NUM_JOBS = int(args.num_jobs)

## Configuration settings
MAX_FILE = -1    ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event




offset = 0

MAX_dR = .1 #Proximity for track matching
s6_chambers = [33, 34, 35, 36, 1, 2] #projections from chambers to sector 6 require some fudging since they aren't contiguous
evt_tree  = TChain('EMTFNtuple/tree') #event tree for analyzing event data


if args.num_jobs:
  try: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/rate_calculation_unpacked%d.root" % (INDEX), 'create')
  except: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/rate_calculation_unpacked%d.root" % (INDEX), 'recreate')
else:
  out_file = TFile('plots/rates/rate_unpacked.root', 'recreate')


#Specify which datasets to use
ZERO_BIAS = True

if ZERO_BIAS:
  base_dirs = ["/eos/user/n/nhurley/EphemeralZeroBias%d/EMTFNtuple_Run3_EphemeralZeroBias%d_data_13p6TeV_Run3Alignment_2022D_v2/" % (i, i) for i in range(21)]
else:
  base_dirs = ["/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run3Alignment_2022C_v5/"]

rate_v_pos_num_1D = np.zeros(shape = (3), dtype=object) #trigger-rate-numerator vs. eta, indexed by BX -2, -1, 0
rate_v_pos_num = np.zeros(shape = (3), dtype=object) # trigger-rate numerator vs. eta-phi, indexed by BX
eta_v_pt = np.zeros(shape = (3), dtype=object) #track occupancy at different pt's and eta's, indexed by BX
eta_v_mode = np.zeros(shape = (3), dtype=object) #track occupancy at different mode's and eta's, indexed by BX
pt_1D = np.zeros(shape = (3), dtype=object) #PT distribution of tracks, indexed by BX, cut by mode (11, 13, 14, 15)
mode_1D = np.zeros(shape = (3), dtype=object) #Mode distribution of tracks, indexed by BX
high_eta_mode = np.zeros(shape = (3), dtype=object) #Mode distribution of tracks only in high-eta region, indexed by BX
high_eta_pt = np.zeros(shape = (3), dtype=object) #PT distribution of tracks only in high-eta region, indexed by BX
high_eta_pos = np.zeros(shape = (3), dtype=object) #occupancy of high-eta tracks in space (this one doesn't really make much sense, redundant with rate_v_pos_num)
rate_v_pos_den = TH2D('rate_v_pos_den', '', 50, -2.4, 2.4, 30, -3.14159, 3.14159) #denominator of all tracks at a given eta-phi
pt_den = np.zeros(shape = (3), dtype=object) #PT distribution of tracks, uncut by mode
eta_rpc = np.zeros(shape = (3), dtype=object)#plot RPC hits in a track vs. eta


#Index by BX -2, -1, 0
mode_den = np.zeros(shape=(2), dtype = object)
for offset in range(-2, 1):
  rate_v_pos_num_1D[offset + 2] = TH1D('rate_v_pos_num_1D_%d' % offset, '', 50, -2.4, 2.4)
  rate_v_pos_num[offset + 2] = TH2D('rate_v_pos_num_%d' % offset, '', 50, -2.4, 2.4, 30, -3.14159, 3.14159)
  eta_v_pt[offset + 2] = TH2D('eta_v_pt_%d' % offset, '', 50, -2.4, 2.4, 70, 0, 70)
  pt_1D[offset + 2] = TH1D('pt_1D_%d' % offset, '', 70, 0, 70)
  mode_1D[offset + 2] = TH1D('mode_1D_%d' % offset, '', 16, 0, 16)
  high_eta_mode[offset + 2] = TH1D('high_eta_mode_%d' % offset, '', 16, 0, 16) 
  high_eta_pt[offset + 2] = TH1D('high_eta_pt_%d' % offset, '', 20, 0, 20) 
  high_eta_pos[offset + 2] = TH2D('high_eta_pos_%d' % offset, '', 30, 2.5, 2.6, 10, -3.14159, 3.14159)
  pt_den[offset + 2] = TH1D('pt_den_%d' % offset, '', 20, 0, 20)
  mode_den[offset + 2] = TH1D('mode_den_%d' % offset, '', 16, 0, 16)
  eta_rpc[offset + 2] = TH1D('eta_rpc_%d' % offset, '', 49, -2.4, 2.4)


bx_den = TH1D('bx_den', '', 3, -2, 1) #get all events in a BX to normalize across BX's in plots listed below
bx_v_pt = TH2D('bx_v_pt', '', 3, -2, 1, 50, 0, 70) #look at pt distribution at different BX's
bx_v_mode = TH2D('bx_v_mode', '', 3, -2, 1, 16, 0, 16) # look at mode distribution at different BX's
high_eta_bx = TH1D('high_eta_bx', '', 3, -2, 1)  #count of high-eta tracks at different bx's
persistance_v_eta = TH1D('persistance_v_eta', '', 49, -2.4, 2.4) #count number of RPC hits that persist between BX-1 and BX0, by eta
bx_rpc = TH1D('bx_rpc', '', 3, -2, 1) #look at distribution of RPC hits in a track vs. BX
dr_rpc = TH1D('dr_rpc', '', 50, 0, .5) #look at distribution between matched RPC hits
persistance_v_eta_rpc = TH1D('persistenc_v_eta_rpc', '', 75, -2.4, 2.4)
persistance_v_eta_csc = TH1D('persistence_v_eta_csc', '', 75, -2.4, 2.4)

matched_all_chambers = TH2D('matched_all_chambers', '', 42, 1, 43, 2 * (1 + 1 + 2 + 2), 0, 12 )


#Set the labels for the rings
for i in range(12):
  ring = '2'
  station = '1'
  if i > 5:
    endcap = '+'
    abs_i = i % 6
  else:
    abs_i = 5 - i
    endcap = '-'
  if abs_i > 0: station = '2'
  if abs_i > 1: station = '3'
  if abs_i > 3: station = '4'
  if abs_i == 3 or abs_i == 5: ring = '3'


  matched_all_chambers.GetYaxis().SetBinLabel(i + 1, 'RE%s%s/%s' % (endcap, station, ring))

#Set the labels for the chambers (accounting for neighbors)
count = 0
xbin = 1
while xbin < 43:
  matched_all_chambers.GetXaxis().SetBinLabel(xbin, str(xbin - count))
  if xbin == 6 or xbin == 13 or xbin == 20 or xbin == 27 or xbin == 34 or xbin == 41:
    xbin += 1
    count += 1
    matched_all_chambers.GetXaxis().SetBinLabel(xbin, "N")
  xbin += 1

matched_chambers = np.zeros(shape=(2, 4, 3), dtype=object)
for endcap in range(2):
  for station in range(4):
    for ring in range(3):
      matched_chambers[endcap][station][ring] = TH1D('matched_chambers_e%d_s%d_r%d' % (endcap, station + 1, ring + 1), '', 36, 1, 37)

#Get all of the event files in the directory
nFiles = 0
#for file in Popen(['ls', folder], stdout=PIPE).communicate()[0].split():


#Get all of the event files in the directory
nFiles = 0


#Flag for breaking loop if we hit max file limit
break_loop = False

#recursivelh access different subdirectories of given folder from above
for base_dir in base_dirs:
  for dirname, dirs, files in os.walk(base_dir):
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

#Get a general count of tracks
denominator = 0
numerator = 0
total_tracks = 0

#Go event by event
for event in range(evt_tree.GetEntries()):
    event_filled = False
    evt_tree.GetEntry(event)

    if event % PRT_EVT == 0:
        if args.index:
            print('rate_calculation_unpacked.py: Processing Job #%d, Event #%d' % (INDEX, event))
        else: print('rate_calculation_unpacked.py: Processing Event #%d' % (event))
    denominator += 1
    run = evt_tree.eventInfo_run[0]

    for track in range(evt_tree.emtfUnpTrack_size):
        track_eta = evt_tree.emtfUnpTrack_eta[track]
        track_phi = evt_tree.emtfUnpTrack_phi[track]
        track_phi *= 3.14159 / 180.0
        if (track_phi > 3.14159): track_phi -= 2*3.14159
        elif (track_phi < -3.14159): track_phi += 2*3.14159

        track_mode = evt_tree.emtfUnpTrack_mode[track]
        track_pt = evt_tree.emtfUnpTrack_pt[track]
        rate_v_pos_den.Fill(track_eta, track_phi)

        total_tracks += 1
        for offset in range(-2, 1):

          #Need to account for the BX shift configuration from August 4th 2022
          if run >= 356798 and evt_tree.emtfUnpTrack_bx[track] != -1 + offset: continue
          elif run < 356798 and evt_tree.emtfUnpTrack_bx[track] != 0 + offset: continue
          bx_den.Fill(offset)
          pt_den[offset + 2].Fill(track_pt)
          if track_pt > 22 - .01 and track_mode in [11, 13, 14, 15]: 
            numerator += 1
            #Add triggered track to rate
            rate_v_pos_num[offset + 2].Fill(track_eta, track_phi)
            rate_v_pos_num_1D[offset + 2].Fill(track_eta)
          if track_mode in [11, 13, 14, 15]:
            #Get pt ditribution of high-quality tracks
            eta_v_pt[offset + 2].Fill(track_eta, track_pt)
            pt_1D[offset + 2].Fill(track_pt)

          #Get mode distribution of all tracks
          mode_1D[offset + 2].Fill(track_mode)
          
          #Get PT/Mode distribution of all tracks w.r.t. PT
          bx_v_pt.Fill(offset, track_pt)
          bx_v_mode.Fill(offset, track_mode)

          #Get rate of tracks with out-of-range eta-values
          if abs(track_eta) > 2.5:
            high_eta_bx.Fill(offset)
            high_eta_mode[offset + 2].Fill(track_mode)
            high_eta_pt[offset + 2].Fill(track_pt)
            high_eta_pos[offset + 2].Fill(track_eta, track_phi)


          for n in range(4):
            hit_ref = eval("evt_tree.emtfUnpTrack_hitref%d[track]" % (n + 1))
            
            if hit_ref < 0: continue

            #Get BX distribution and eta of RPC hits in tracks
            if evt_tree.emtfUnpHit_type[hit_ref] == 2: 
              bx_rpc.Fill(offset)
              eta_rpc[offset + 2].Fill(track_eta)

    for hit_ref in range(evt_tree.emtfUnpHit_size):

      #Need to account for the BX shift configuration from August 4th 2022
      #Get hit from BX-1
      if run >= 356798 and evt_tree.emtfUnpHit_bx[hit_ref] != -2: continue
      elif run < 356798 and evt_tree.emtfUnpHit_bx[hit_ref] != -1: continue
      #Get hit information
      hit_theta = evt_tree.emtfUnpHit_emtf_theta[hit_ref]
      if hit_theta < 0: continue
      hit_theta = hit_theta * (45. - 8.5) / 128. + 8.5
      hit_theta *= 3.14159 / 180.0
      hit_eta = -math.log(math.tan(hit_theta/2.0))

      if evt_tree.emtfUnpHit_endcap[hit_ref] == -1: hit_eta *= -1
      hit_sector = evt_tree.emtfUnpHit_sector[hit_ref]

      #Convert hit phi to global coordinates
      hit_phi_loc = evt_tree.emtfUnpHit_emtf_phi[hit_ref] / 60.0 - 22.
      hit_phi_glob = hit_phi_loc + 15. + (60. * (hit_sector - 1))
      if hit_phi_glob > 180.0: hit_phi_glob -= 360.
      hit_phi = hit_phi_glob * 3.14159/180.

      #Only take RPC hits
      if evt_tree.emtfUnpHit_type[hit_ref] != 2: continue
      #Now we match RPC hit in BX-1 to a hit in BX0
      for hit_ref2 in range(evt_tree.emtfUnpHit_size):

        #Need to account for the BX shift configuration from August 4th 2022
        #Get hit from BX0
        if run >= 356798 and evt_tree.emtfUnpHit_bx[hit_ref2] != -1: continue
        elif run < 356798 and evt_tree.emtfUnpHit_bx[hit_ref2] != 0: continue

        #Get hit coordinates
        hit_theta2 = evt_tree.emtfUnpHit_emtf_theta[hit_ref2]
        if hit_theta2 < 0: continue
        hit_theta2 = hit_theta2 * (45. - 8.5) / 128. + 8.5
        hit_theta2 *= 3.14159 / 180.0
        hit_eta2 = -math.log(math.tan(hit_theta2/2.0))
        if evt_tree.emtfUnpHit_endcap[hit_ref2] == -1: hit_eta2 *= -1
        hit_sector2 = evt_tree.emtfUnpHit_sector[hit_ref2]

        #Convert hit phi to global coordinates
        hit_phi_loc2 = evt_tree.emtfUnpHit_emtf_phi[hit_ref2] / 60.0 - 22.
        hit_phi_glob2 = hit_phi_loc2 + 15. + (60. * (hit_sector2 - 1))
        if hit_phi_glob2 > 180.0: hit_phi_glob2 -= 360.
        hit_phi2 = hit_phi_glob2 * 3.14159/180.

        #only take RPC hits
        if evt_tree.emtfUnpHit_type[hit_ref2] != 2: continue

        #Matched two hits by distance
        dr = h.CalcDR(hit_eta, hit_phi, hit_eta2, hit_phi2)
        #Course matching of hit distance to see distribution
        if dr > 2 * MAX_dR: continue
        dr_rpc.Fill(dr)

        #Match hits more strictly
        if dr <= .05:

          #Match hits by chamber, ensure type RPC
          if evt_tree.emtfUnpHit_type[hit_ref2] == 2 and evt_tree.emtfUnpHit_type[hit_ref] == 2:
            if evt_tree.emtfUnpHit_chamber[hit_ref2] == evt_tree.emtfUnpHit_chamber[hit_ref]:
              endcap = evt_tree.emtfUnpHit_endcap[hit_ref]
              ring = (evt_tree.emtfUnpHit_ring[hit_ref] - 1) % 3
              station = evt_tree.emtfUnpHit_station[hit_ref] - 1
              chamber = evt_tree.emtfUnpHit_chamber[hit_ref]
              
              if not evt_tree.emtfUnpHit_neighbor[hit_ref]:
                persistance_v_eta_rpc.Fill(hit_eta)
                matched_chambers[endcap == 1][station][ring].Fill(chamber)


              #Get station-ring y-coordinate
              ring += 1
              station += 1
              index = station - 1 

              if ring == 3: index += 1
              if station == 4: index += 1

              #Endcap 1 is 6-11
              #Endcap -1 is 5-0 (reversed)
              if endcap == 1: index += 6
              elif endcap == -1: index = 5 - index


              #Fill chamber or neighbor chamber
              if not evt_tree.emtfUnpHit_neighbor[hit_ref]:
                chamber_bin = chamber + ((chamber - 1) / 6)
                matched_all_chambers.Fill(chamber_bin, index)
              else:
                neighbor_bin = (evt_tree.emtfUnpHit_sector[hit_ref]) * 7
                matched_all_chambers.Fill(neighbor_bin, index)

          #Match csc hits in BX to compare to RPC hits
          elif evt_tree.emtfUnpHit_type[hit_ref2] == 1 and evt_tree.emtfUnpHit_type[hit_ref] == 1:
            if evt_tree.emtfUnpHit_chamber[hit_ref2] == evt_tree.emtfUnpHit_chamber[hit_ref]:
              persistance_v_eta_csc.Fill(hit_eta)
            

print("Total tracks: " + str(total_tracks))
print("Tagged events: " + str(numerator))
print("Total events: "+ str(denominator))

#Write all of our plots
rate_v_pos_den.Write()
for offset in range(-2, 1):
  rate_v_pos_num_1D[offset + 2].Write()
  rate_v_pos_num[offset + 2].Write()
  eta_v_pt[offset + 2].Write()
  pt_1D[offset + 2].Write()
  mode_1D[offset + 2].Write()
  high_eta_mode[offset + 2].Write()
  high_eta_pt[offset + 2].Write()
  high_eta_pos[offset + 2].Write()
  pt_den[offset + 2].Write()
  mode_den[offset + 2].Write()
  eta_rpc[offset + 2].Write()

high_eta_bx.Write()
bx_den.Write()
bx_rpc.Write()
persistance_v_eta_rpc.Write()
persistance_v_eta_csc.Write()
dr_rpc.Write()
matched_all_chambers.Write()
for endcap in range(2):
  for station in range(4):
    for ring in range(3):
      matched_chambers[endcap][station][ring].Write()


for bx in range(-2, 1):
  for pt in range(0, 70):
    bx_v_pt.SetBinContent(bx + 3, pt + 1, bx_v_pt.GetBinContent(bx + 3, pt + 1) / bx_den.GetBinContent(bx + 3))
bx_v_pt.Write()

for bx in range(-2, 1):
  for mode in range(0, 16):
    bx_v_mode.SetBinContent(bx + 3, mode + 1, bx_v_mode.GetBinContent(bx + 3, mode + 1) / bx_den.GetBinContent(bx + 3))
bx_v_mode.Write()


# canvas = TCanvas(rate_v_pos_num.GetName() , rate_v_pos_num.GetName(), 700,700)

# gStyle.SetOptStat(0)

# rate_v_pos_num.Draw("colz")
# rate_v_pos_num.GetXaxis().SetTitle("#eta")
# rate_v_pos_num.GetYaxis().SetTitle("#phi")
# cms_label = cms_latex()
# header = head()

# canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/rates/" + rate_v_pos_num.GetName() + "%d.pdf" % (offset))


# rate_v_pos_num.Write()
# rate_v_pos_den.Write()
# rate_v_pos_num.Divide(rate_v_pos_den)
# rate_v_pos_num.SetName('rate_v_pos_tot')
# rate_v_pos_num.Write()

# rate_v_pos_num_1D.Write()

# rate = float(numerator) / float(denominator)
# print("Rate is " + "%1.8f" % (rate))


del out_file
