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


''' This script exists to analyze the rates of EMTF tracks'''

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
evt_tree  = TChain('EMTFNtuple/tree') #event tree for analyzing event data


if args.num_jobs:
  try: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/rate_calculation%d.root" % (INDEX), 'create')
  except: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/rate_calculation%d.root" % (INDEX), 'recreate')
else:
  out_file = TFile('plots/rate.root', 'recreate')


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
  eta_rpc[offset + 2] = TH1D('eta_rpc_%d' % offset, '', 49, -2.4, 2.4)


bx_den = TH1D('bx_den', '', 3, -2, 1) #get all events in a BX to normalize across BX's in plots listed below
bx_v_pt = TH2D('bx_v_pt', '', 3, -2, 1, 50, 0, 70) #look at pt distribution at different BX's
bx_v_mode = TH2D('bx_v_mode', '', 3, -2, 1, 16, 0, 16) # look at mode distribution at different BX's
high_eta_bx = TH1D('high_eta_bx', '', 3, -2, 1)  #count of high-eta tracks at different bx's
persistance_v_eta = TH1D('persistance_v_eta', '', 49, -2.4, 2.4) #count number of RPC hits that persist between BX-1 and BX0, by eta
bx_rpc = TH1D('bx_rpc', '', 3, -2, 1) #look at distribution of RPC hits in a track vs. BX



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
            print('rate_calculation.py: Processing Job #%d, Event #%d' % (INDEX, event))
        else: print('rate_calculation.py: Processing Event #%d' % (event))
        
    #Count every event
    denominator += 1
    run = evt_tree.eventInfo_run[0]

    for track in range(evt_tree.emtfTrack_size):
        track_eta = evt_tree.emtfTrack_eta[track]
        track_phi = evt_tree.emtfTrack_phi[track]
        track_phi *= 3.14159 / 180.0
        if (track_phi > 3.14159): track_phi -= 2*3.14159
        elif (track_phi < -3.14159): track_phi += 2*3.14159

        track_mode = evt_tree.emtfTrack_mode[track]
        track_pt = evt_tree.emtfTrack_pt[track]
        rate_v_pos_den.Fill(track_eta, track_phi)

        
        total_tracks += 1
        for offset in range(-2, 1):

          #Need to account for the BX shift configuration from August 4th 2022
          if run >= 356798 and evt_tree.emtfTrack_bx[track] != -1 + offset: continue
          elif run < 356798 and evt_tree.emtfTrack_bx[track] != 0 + offset: continue
          
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
            hit_ref = eval("evt_tree.emtfTrack_hitref%d[track]" % (n + 1))
            if hit_ref < 0: continue
            
            #Get BX distribution and eta of RPC hits in tracks
            if evt_tree.emtfHit_type[hit_ref] == 2: 
              bx_rpc.Fill(offset)
              eta_rpc[offset + 2].Fill(track_eta)

    for hit_ref in range(evt_tree.emtfHit_size):

      #Need to account for the BX shift configuration from August 4th 2022
      #Get hit from BX-1
      if run >= 356798 and evt_tree.emtfHit_bx[hit_ref] != -2: continue
      elif run < 356798 and evt_tree.emtfHit_bx[hit_ref] != -1: continue


      #Get hit information
      hit_theta = evt_tree.emtfHit_emtf_theta[hit_ref]
      if hit_theta < 0: continue
      hit_theta = hit_theta * (45. - 8.5) / 128. + 8.5
      hit_theta *= 3.14159 / 180.0
      hit_eta = -math.log(math.tan(hit_theta/2.0))

      if evt_tree.emtfHit_endcap[hit_ref] == -1: hit_eta *= -1
      hit_sector = evt_tree.emtfHit_sector[hit_ref]

      #Convert hit phi to global coordinates
      hit_phi_loc = evt_tree.emtfHit_emtf_phi[hit_ref] / 60.0 - 22.
      hit_phi_glob = hit_phi_loc + 15. + (60. * (hit_sector - 1))
      if hit_phi_glob > 180.0: hit_phi_glob -= 360.
      hit_phi = hit_phi_glob * 3.14159/180.

      #Only take RPC hits
      if evt_tree.emtfHit_type[hit_ref] != 2: continue


      #Now we match RPC hit in BX-1 to a hit in BX0
      for hit_ref2 in range(evt_tree.emtfHit_size):

        #Need to account for the BX shift configuration from August 4th 2022
        #Get hit from BX0
        if run >= 356798 and evt_tree.emtfHit_bx[hit_ref2] != -1: continue
        elif run < 356798 and evt_tree.emtfHit_bx[hit_ref2] != 0: continue

        #Get hit coordinates
        hit_theta2 = evt_tree.emtfHit_emtf_theta[hit_ref2]
        if hit_theta2 < 0: continue
        hit_theta2 = hit_theta2 * (45. - 8.5) / 128. + 8.5
        hit_theta2 *= 3.14159 / 180.0
        hit_eta2 = -math.log(math.tan(hit_theta2/2.0))
        if evt_tree.emtfHit_endcap[hit_ref2] == -1: hit_eta2 *= -1
        hit_sector2 = evt_tree.emtfHit_sector[hit_ref2]

        #Convert hit phi to global coordinates
        hit_phi_loc2 = evt_tree.emtfHit_emtf_phi[hit_ref2] / 60.0 - 22.
        hit_phi_glob2 = hit_phi_loc2 + 15. + (60. * (hit_sector2 - 1))
        if hit_phi_glob2 > 180.0: hit_phi_glob2 -= 360.
        hit_phi2 = hit_phi_glob2 * 3.14159/180.

        #only take RPC hits
        if evt_tree.emtfHit_type[hit_ref2] != 2: continue


        #Matched two hits by distance
        dr = h.CalcDR(hit_eta, hit_phi, hit_eta2, hit_phi2)
        if dr > 2 * MAX_dR: continue

        persistance_v_eta.Fill(hit_eta)     

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
  eta_rpc[offset + 2].Write()

high_eta_bx.Write()
bx_den.Write()
bx_rpc.Write()
persistance_v_eta.Write()


for bx in range(-2, 1):
  for pt in range(0, 70):
    bx_v_pt.SetBinContent(bx + 3, pt + 1, bx_v_pt.GetBinContent(bx + 3, pt + 1) / bx_den.GetBinContent(bx + 3))
bx_v_pt.Write()

for bx in range(-2, 1):
  for mode in range(0, 16):
    bx_v_mode.SetBinContent(bx + 3, mode + 1, bx_v_mode.GetBinContent(bx + 3, mode + 1) / bx_den.GetBinContent(bx + 3))
bx_v_mode.Write()
  
del out_file
