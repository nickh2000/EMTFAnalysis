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
MAX_FILE =  -1         ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1
       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event

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

#plot a track's delta-phi between stations as a function of phi
#encap, 2-station-tuple, sector (0 index is all sectors)
d_phi_plots = np.zeros(shape = (2, 6, 7, 3), dtype=object)

#do same but with chambers instead of sectors
d_phi_chamber = np.zeros(shape = (2, 3, 6, 36, 3), dtype=object)
num_bins = 75

for endcap in range(2):
    for d in range(6):
    #index 0 includes all sectors to get a global dphi offset
        d_phi_plots[endcap][d][0][0] = TH1D('d_phi_EMTF_%s_%d_all' % (str(d_phi_index_map[d]), endcap),  '', 128, -64, 64)

        #then plot by sector
        for sector in range(6):
            d_phi_plots[endcap][d][sector + 1][0] = TH1D('d_phi_EMTF_%s_%d_%d' % (str(d_phi_index_map[d]), endcap, sector + 1),  '', 128, -64, 64)
            d_phi_plots[endcap][d][sector + 1][0].GetYaxis().SetTitle('Instances')
            d_phi_plots[endcap][d][sector + 1][0].GetXaxis().SetTitle("#Delta#phi_{%d%d}" % (d_phi_index_map[d][0], d_phi_index_map[d][1]))

            d_phi_plots[endcap][d][sector + 1][0].SetTitle("#Delta#phi_{%d%d}, Endcap %d, Sector %d" %
                (d_phi_index_map[d][0], 
                d_phi_index_map[d][1],
                1 if endcap else -1, 
                sector + 1))

        #also do plots for only negative/positive muons
        for sector in range(7):
            for charge in range(1,3):
                d_phi_plots[endcap][d][sector][charge] = TH1D('d_phi_EMTF_%s_%d_%d_c%d' % (str(d_phi_index_map[d]), endcap, sector, charge - 1),  '', 128, -64, 64)

        #do same for chambers instead of sectors
        for chamber in range(36):
            for ring in range(3):
                for charge in range(3):
                    d_phi_chamber[endcap][ring][d][chamber][charge] = TH1D('d_phi_chamber_%s_%d_%d_%d_c%d' % (str(d_phi_index_map[d]), endcap, ring, chamber, charge),  '', 128, -64, 64)

#Just save arguments as variables 
if args.num_jobs:
    INDEX = int(args.index)
    NUM_JOBS = int(args.num_jobs)


#define output files
if args.num_jobs:
  try: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/station_station_phi_custom_chamber%d.root" % (INDEX), 'create')
  except: out_file =  TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/station_station_phi_custom_chamber%d.root" % (INDEX), 'recreate')
else:
  out_file = TFile('plots/high_eta_region/high_eta_study_Custom.root', 'recreate')


#Determine whether to use Ideal, Run3, Run2, or customized-by-me geometry look-up tables in our EMTFNTuple Re-Emulation
#These folders are collections of EMTFNTuples, with cscSegments, re-emulated&unpacked hits, tracks enabled, ,and offline-reconstructed muons
#CMSSW Produceder for EMTFNTuples is found here: https://github.com/eyigitba/EMTFTools/tree/master/EMTFNtuple


IDEAL = False
RUN3 = False
CUSTOM = True

if IDEAL:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_SingleMuon_data_13p6TeV_idealAlignment/220902_142649/0000/"
elif RUN3:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run3Alignment_2022C_v5/220925_185603/0000/"
elif CUSTOM:
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_phlutv3data_2022C_v1/221020_181346/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v9/221020_120011/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v8/221019_170559/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_init_v1/221021_100244/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_init_v2/221022_102914/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_reverse_v2/221023_205023/"
    # folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_reverse_v3/221026_134022/0000/"
    # folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_reverse_v4/221026_174204/0000/"
    #folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v11/221028_095517/0000/"
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v10/221027_083620/0000/"
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_chamber_v1/221109_095923/0000/"
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_chamber_v3/221109_220323/0000/"
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_chamber_v4/221110_130238/0000/"
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_CustomAlignment_2022C_v15/221111_132410/0000/"
else:
    folder = "/eos/user/n/nhurley/Muon/EMTFNtuple_Run3_Muon_data_13p6TeV_Run2Alignment_2022C_v2/220920_155151/0000/"
    
#Get all of the event files in the directory
nFiles = 0


#recursively access different subdirectories of given folder from above
break_loop = False
for dirname, dirs, files in os.walk(folder):
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
        print ('* Loading file #%s: %s' % (nFiles, file))
        evt_tree.Add(file_name)
        if nFiles == MAX_FILE: break_loop = True


#Iterate through events in EMTFNtuple
for event in range(evt_tree.GetEntries()):

    evt_tree.GetEntry(event)
    if event == MAX_EVT: break
    if event % PRT_EVT == 0:
        if args.index:
            print('station_station_phi.py: Processing Job #%d, Event #%d' % (INDEX, event))
        else: print('station_station_phi.py: Processing Event #%d' % (event))

    for track in range(evt_tree.emtfTrack_size):
    # #Station-Station DPhi plots
        # #Find which station-station transsitions exist based on the mode
        mode = evt_tree.emtfTrack_mode[track]
        if not mode in d_phi_mode_map: continue

        #list of different transitions by mode
        for id in d_phi_mode_map[mode]:
            

            #Get the two station hits
            first_station = d_phi_index_map[id][0]
            second_station = d_phi_index_map[id][1]

            #Get dphi between the two hits at the different stations
            hitref_1 = eval("evt_tree.emtfTrack_hitref%d[track]" % (first_station))
            hitref_2 = eval("evt_tree.emtfTrack_hitref%d[track]" % (second_station))

            d_phi = evt_tree.emtfHit_emtf_phi[hitref_2] - evt_tree.emtfHit_emtf_phi[hitref_1]

            #Fill d_phi by phi for all sectors, both charges
            d_phi_plots[evt_tree.emtfTrack_endcap[track] == 1][id][0][0].Fill(d_phi)
            

            #Fill dphi for specific sector, all charges
            d_phi_plots[evt_tree.emtfTrack_endcap[track] == 1][id][evt_tree.emtfTrack_sector[track]][0].Fill(d_phi)

            if evt_tree.emtfHit_ring[hitref_1] == 4:
                ring = 1
            else: ring = evt_tree.emtfHit_ring[hitref_1]

            #Fill dphi for all chambers, all charges
            d_phi_chamber[evt_tree.emtfTrack_endcap[track] == 1][ring-1][id][evt_tree.emtfHit_chamber[hitref_1] - 1][0].Fill(d_phi)
            charge = (evt_tree.emtfTrack_q[track] == 1) + 1

            #Fill dphi for specific charge, specific sector
            d_phi_plots[evt_tree.emtfTrack_endcap[track] == 1][id][evt_tree.emtfTrack_sector[track]][charge].Fill(d_phi)
            
            #Fill dphi for specific charge, specific sector
            d_phi_chamber[evt_tree.emtfTrack_endcap[track] == 1][ring - 1][id][evt_tree.emtfHit_chamber[hitref_1] - 1][charge].Fill(d_phi)


#Write station-station phi plots to output file
for endcap in range(2):
    for d in range(6):
        for sector in range(7):
            for charge in range(3):
                d_phi_plots[endcap][d][sector][charge].Write()
        for chamber in range(36):
            for ring in range(3):
                for charge in range(3):
                    d_phi_chamber[endcap][ring][d][chamber][charge].Write()


del out_file