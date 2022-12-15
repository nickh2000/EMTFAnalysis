from cmath import tan
import sys
import math
from subprocess import Popen,PIPE
from turtle import clear
from ROOT import *
import numpy as np
from array import *
import Helper as h
import argparse
import os


'''This script exists to find the distrubutions of delta-phi's between stations at each sector
    Script is used to determine extent of misalignment
    output files are used to corect misalignment in shift_alignment_from_data.py (set station_station_phi_data directory string in file)
'''


#Configure input arguments
#num_jobs: split our input files into n different processes
#index: the index of the current running process, defines which chunk of input files to fetch
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--num_jobs", required=False)
parser.add_argument("-i", "--index", required = False)
args = parser.parse_args()

## Configuration settings
MAX_FILE =  -1        ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1
       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event

## HISTOGRAMS 


#Just save arguments as variables 
if args.num_jobs:
    INDEX = int(args.index)
    NUM_JOBS = int(args.num_jobs)


#define output files
if args.num_jobs:
  try: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/2018compare%d.root" % (INDEX), 'create')
  except: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/2018compare%d.root" % (INDEX), 'recreate')
else:
  out_file = TFile('plots/2018compare/2018compare.root', 'recreate')


#Determine whether to use Ideal, Run3, Run2, or customized-by-me geometry look-up tables in our EMTFNTuple Re-Emulation
#These folders are collections of EMTFNTuples, with cscSegments, re-emulated&unpacked hits, tracks enabled, ,and offline-reconstructed muons
#CMSSW Produceder for EMTFNTuples is found here: https://github.com/eyigitba/EMTFTools/tree/master/EMTFNtuple
scale_pt_temp = [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 45, 60, 75, 100, 140, 160, 180, 200, 250, 300, 500, 1000] #high-pt range
scale_pt  = array('d', scale_pt_temp)

scale_pt_temp_2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 55, 60] #low-pt range
scale_pt_2  = array('d', scale_pt_temp_2)

pt_2018 = TH1D('pt_EMTF_2018',  '', len(scale_pt_temp) - 1,  scale_pt)
pt_2022 = TH1D('pt_EMTF_2022',  '', len(scale_pt_temp) - 1,  scale_pt)

pt_2_2018 = np.zeros(shape=(3), dtype=object)
pt_2_2022 = np.zeros(shape=(3), dtype=object)
for q in range(3):
    pt_2_2018[q] = TH1D('pt_2_%d_EMTF_2018' % q,  '', len(scale_pt_temp_2) - 1,  scale_pt_2)
    pt_2_2022[q] = TH1D('pt_2_%d_EMTF_2022' % q,  '', len(scale_pt_temp_2) - 1,  scale_pt_2)

eta_2018 = TH1D('eta_EMTF_2018',  '', 50,  -2.5, 2.5)
eta_2022 = TH1D('eta_EMTF_2022',  '', 50,  -2.5, 2.5)

phi_2018 = TH1D('phi_EMTF_2018',  '', 50,  -3.14159, 3.14159)
phi_2022 = TH1D('phi_EMTF_2022',  '', 50,  -3.14159, 3.14159)



pt_eta_2018 = TH2D('pt-eta_EMTF_2018', '', 50, -2.5, 2.5, len(scale_pt_temp_2) - 1, scale_pt_2)
pt_eta_2022 = TH2D('pt-eta_EMTF_2022', '', 50, -2.5, 2.5, len(scale_pt_temp_2) - 1, scale_pt_2)

pt_phi_2018 = TH2D('pt-phi_EMTF_2018', '', 50, -3.14159, 3.14159, len(scale_pt_temp_2) - 1, scale_pt_2)
pt_phi_2022 = TH2D('pt-phi_EMTF_2022', '', 50, -3.14159, 3.14159, len(scale_pt_temp_2) - 1, scale_pt_2)

rate_eta_2018 = TH1D('rate-eta_EMTF_2018', '', 50, -2.5, 2.5)
rate_eta_2022 = TH1D('rate-eta_EMTF_2022', '', 50, -2.5, 2.5)

rate_phi_2018 = TH1D('rate-phi_EMTF_2018', '', 50, -3.14159, 3.14159)
rate_phi_2022 = TH1D('rate-phi_EMTF_2022', '', 50, -3.14159, 3.14159)

rate_eta_phi_2018 = TH2D('rate-eta-phi_EMTF_2018', '', 50, -2.5, 2.5, 50, -3.14159, 3.14159)
rate_eta_phi_2022 = TH2D('rate-eta-phi_EMTF_2022', '', 50, -2.5, 2.5, 50, -3.14159, 3.14159)

pt_eta_bin_2018 = np.zeros(shape=(3, 3), dtype=object)
pt_eta_bin_2022 = np.zeros(shape=(3, 3), dtype=object)

for i in range(3):
    for q in range(3):
        pt_eta_bin_2018[i][q] = TH1D('pt_2_%d_%d_EMTF_2018' % (i, q),  '', len(scale_pt_temp_2) - 1,  scale_pt_2)
        pt_eta_bin_2022[i][q] = TH1D('pt_2_%d_%d_EMTF_2022' % (i, q),  '', len(scale_pt_temp_2) - 1,  scale_pt_2)


IN_2018 = ["/eos/user/n/nhurley/EphemeralZeroBias%d/EMTFNtuple_Collisions18_EphemeralZeroBias%d_data_13TeV_2018D_v1/" % (d, d) for d in range(20)]
IN_2022 = ["/eos/user/n/nhurley/EphemeralZeroBias%d/EMTFNtuple_Collisions22_EphemeralZeroBias%d_data_13TeV_2022F_v1/" % (d, d) for d in range(20)]
    
#Get all of the event files in the directory
nFiles = 0


evt_tree_2018  = TChain('EMTFNtuple/tree')
evt_tree_2022  = TChain('EMTFNtuple/tree')

#recursively access different subdirectories of given folder from above
break_loop = False
for dir in IN_2018:
    for dirname, dirs, files in os.walk(dir):
        if break_loop: break

        if args.num_jobs and args.index:
            file_list = files[INDEX * len(files) / NUM_JOBS : (INDEX + 1) * len(files) / NUM_JOBS]
        else:
            file_list = files
        for file in file_list:
            if break_loop: break
            if not file.find('.root') == len(file) - 5: continue
            file_name = "%s/%s" % (dirname, file)
            nFiles   += 1
            print ('* Loading file #%s: %s' % (nFiles, file_name))
            evt_tree_2018.Add(file_name)
            if nFiles == MAX_FILE: break_loop = True



nFiles = 0
break_loop = False
for dir in IN_2022:
    for dirname, dirs, files in os.walk(dir):
        if break_loop: break

        if args.num_jobs and args.index:
            file_list = files[INDEX * len(files[0:MAX_FILE]) / NUM_JOBS : (INDEX + 1) * len(files[:MAX_FILE]) / NUM_JOBS]
        else:
            file_list = files[0:MAX_FILE]

        for file in file_list:
            if break_loop: break
            if not file.find('.root') == len(file) - 5: continue
            file_name = "%s/%s" % (dirname, file)
            nFiles   += 1
            print ('* Loading file #%s: %s' % (nFiles, file_name))
            evt_tree_2022.Add(file_name)
            if nFiles == MAX_FILE: break_loop = True



for event in range(evt_tree_2018.GetEntries()):
    evt_tree_2018.GetEntry(event)

    if event % PRT_EVT == 0:
        if args.index:
            print('eff_EMTF.py: Processing Job #%d, Year 2018, Event #%d' % (INDEX, event))
        else: print('eff_EMTF.py: Processing Year 2018, Event #%d' % (event))

    if event == MAX_EVT: break
    
    for track in range(evt_tree_2018.emtfUnpTrack_size):
        eta = evt_tree_2018.emtfUnpTrack_eta[track]
        pt = evt_tree_2018.emtfUnpTrack_pt[track]
        mode = evt_tree_2018.emtfUnpTrack_mode[track]

        pt_2018.Fill(evt_tree_2018.emtfUnpTrack_pt[track])
        eta_2018.Fill(evt_tree_2018.emtfUnpTrack_eta[track])
        phi_2018.Fill(evt_tree_2018.emtfUnpTrack_phi[track])
        pt_eta_2018.Fill(evt_tree_2018.emtfUnpTrack_eta[track], evt_tree_2018.emtfUnpTrack_pt[track])
        pt_phi_2018.Fill(evt_tree_2018.emtfUnpTrack_phi[track], evt_tree_2018.emtfUnpTrack_pt[track])


        if mode ==15:
            quality = 15
        elif mode == 14:
            quality = 14
        elif mode == 13:
            quality = 13
        elif mode == 12:
            quality = 7
        elif mode == 11:
            quality = 12
        elif mode == 10:
            quality = 10
        elif mode == 9:
            quality = 9
        elif mode == 7:
            quality = 11
        elif mode == 6:
            quality = 6
        elif mode == 5:
            quality = 5
        elif mode == 3:
            quality = 4
        else:
            quality = 0
    
        qual_index = (quality >= 4) + (quality >= 8) + (quality >= 12)



        for i in range(qual_index):
            pt_2_2018[i].Fill(pt)
            if abs(eta) < 1.6:
                pt_eta_bin_2018[0][i].Fill(pt)
            elif abs(eta) < 2.1:
                pt_eta_bin_2018[1][i].Fill(pt)
            elif abs(eta) < 2.5:
                pt_eta_bin_2018[2][i].Fill(pt)


        if evt_tree_2018.emtfUnpTrack_pt[track] > 22:
            #  and evt_tree_2018.emtfUnpTrack_mode[track] in [11, 13, 14, 15]:
            phi = evt_tree_2018.emtfUnpTrack_phi[track]
            if phi > 3.14159: phi -= 2 * 3.14159
            elif phi < -3.14159: phi += 2 * 3.14159
            rate_eta_2018.Fill(evt_tree_2018.emtfUnpTrack_eta[track])
            rate_phi_2018.Fill(phi)
            rate_eta_phi_2018.Fill(evt_tree_2018.emtfUnpTrack_eta[track], phi)


for event in range(evt_tree_2022.GetEntries()):
    evt_tree_2022.GetEntry(event)

    if event == MAX_EVT: break

    if event % PRT_EVT == 0:
        if args.index:
            print('eff_EMTF.py: Processing Job #%d, Year 2022, Event #%d' % (INDEX, event))
        else: print('eff_EMTF.py: Processing Year 2022, Event #%d' % (event))
    
    for track in range(evt_tree_2022.emtfUnpTrack_size):
        eta = evt_tree_2022.emtfUnpTrack_eta[track]
        pt = evt_tree_2022.emtfUnpTrack_pt[track]
        mode =  evt_tree_2022.emtfUnpTrack_mode[track]

        pt_2022.Fill(evt_tree_2022.emtfUnpTrack_pt[track])
        eta_2022.Fill(evt_tree_2022.emtfUnpTrack_eta[track])
        phi_2022.Fill(evt_tree_2022.emtfUnpTrack_phi[track])
        pt_eta_2022.Fill(evt_tree_2022.emtfUnpTrack_eta[track], evt_tree_2022.emtfUnpTrack_pt[track])
        pt_phi_2022.Fill(evt_tree_2022.emtfUnpTrack_phi[track], evt_tree_2022.emtfUnpTrack_pt[track])
        

        if mode ==15:
            quality = 15
        elif mode == 14:
            quality = 14
        elif mode == 13:
            quality = 13
        elif mode == 12:
            quality = 7
        elif mode == 11:
            quality = 12
        elif mode == 10:
            quality = 10
        elif mode == 9:
            quality = 9
        elif mode == 7:
            quality = 11
        elif mode == 6:
            quality = 6
        elif mode == 5:
            quality = 5
        elif mode == 3:
            quality = 4
        else:
            quality = 0
    
        qual_index = (quality >= 4) + (quality >= 8) + (quality >= 12)

        for i in range(qual_index):
            pt_2_2022[i].Fill(pt)
            if abs(eta) < 1.6:
                pt_eta_bin_2022[0][i].Fill(pt)
            elif abs(eta) < 2.1:
                pt_eta_bin_2022[1][i].Fill(pt)
            elif abs(eta) < 2.5:
                pt_eta_bin_2022[2][i].Fill(pt)

        if evt_tree_2022.emtfUnpTrack_pt[track] > 22:
            #  and evt_tree_2022.emtfUnpTrack_mode[track] in [11, 13, 14, 15]:
            phi = evt_tree_2022.emtfUnpTrack_phi[track]
            if phi > 3.14159: phi -= 2 * 3.14159
            elif phi < -3.14159: phi += 2 * 3.14159
            rate_eta_2022.Fill(evt_tree_2022.emtfUnpTrack_eta[track])
            rate_phi_2022.Fill(phi)
            rate_eta_phi_2022.Fill(evt_tree_2022.emtfUnpTrack_eta[track], phi)

pt_2018.Write()
eta_2018.Write()
phi_2018.Write()
pt_eta_2018.Write()
pt_phi_2018.Write()
rate_eta_2018.Write()
rate_phi_2018.Write()
rate_eta_phi_2018.Write()


pt_2022.Write()
eta_2022.Write()
phi_2022.Write()
pt_eta_2022.Write()
pt_phi_2022.Write()
rate_phi_2022.Write()
rate_eta_2022.Write()
rate_eta_phi_2022.Write()

for q in range(3):
    pt_2_2018[q].Write()
    pt_2_2022[q].Write()


for i in range(3):
    for q in range(3):
        pt_eta_bin_2018[i][q].Write()
        pt_eta_bin_2022[i][q].Write()

del out_file