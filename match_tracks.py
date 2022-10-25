import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
import Helper as h

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

## Configuration settings
USE_EMUL = True    ## Use emulated L1T muons instead of unpacked
MAX_FILE =  -1       ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event

REQ_BX0    = True  ## Require L1T muon to be in BX 0  ## Require a final uGMT candidate, not just a TF muon
## REQ_Z      = False ## Require tag and probe muon to satisfy 81 < mass < 101 GeV (not yet implemented - AWB 25.04.2019)
REQ_TIGHT = False 
MEDIUM_ONLY = False ## Require that muons be of medium quality but not of tight quality
REQ_HLT = False

MAX_dR  = 0.1     ## Maximum dR for L1T-offline matching
TAG_ISO = 0.15   ## Maximum relative isolation for tag muon
TAG_PT  = 26    ## Minimum offline pT for tag muon
PRB_PT  = 26    ## Minimum offline pT for probe muon
TRG_PT  = 22    ## Minimum L1T pT for probe muon


evt_tree  = TChain('EMTFNtuple/tree')
folder = "/eos/cms/store/user/eyigitba/emtf/L1Ntuples/Run3/EfficiencyStudies/357479/EMTFNtuples/Muon/EMTFNtuple_Muon_run3BDT_357479_Run2Alignment_20220825_114518/220825_164535/0000/"
out_file = TFile('plots/L1T_eff_Pt_EMTF.root', 'recreate')

#Get all of the event files in the directory
nFiles = 0
for file in Popen(['ls', folder], stdout=PIPE).communicate()[0].split():
    if not '.root' in file: continue
    file_name = "%s%s" % (folder, file)
    nFiles   += 1
    print ('* Loading file #%s: %s' % (nFiles, file))
    evt_tree.Add(file_name)
    if nFiles == MAX_FILE: break

mismatched_tracks = 0
num_tracks = 0
for event in range(evt_tree.GetEntries()):
    evt_tree.GetEntry(event)
    if event % PRT_EVT == 0:
        print('Processing Event #%d' % (event))

    num_tracks += evt_tree.emtfTrack_size
    for i in range(evt_tree.emtfTrack_size):
        for j in range(evt_tree.emtfUnpTrack_size):
            if evt_tree.emtfTrack_mode[i] == evt_tree.emtfUnpTrack_mode[j]:
                if evt_tree.emtfTrack_phi[i] == evt_tree.emtfUnpTrack_phi[j]:
                    if evt_tree.emtfTrack_sector[i] == evt_tree.emtfUnpTrack_sector[j]:
                        if evt_tree.emtfTrack_bx[i] == evt_tree.emtfUnpTrack_bx[j]:
                            if abs(evt_tree.emtfTrack_eta[i] - evt_tree.emtfUnpTrack_eta[j]) >= .001:
                                mismatched_tracks += 1
                                print("Track 1: Mode=%f Eta=%f Phi=%f Pt=%f" % (evt_tree.emtfTrack_mode[i], evt_tree.emtfTrack_eta[i], evt_tree.emtfTrack_phi[i], evt_tree.emtfTrack_pt[i]))
                                print("Track 2: Mode=%f Eta=%f Phi=%f Pt=%f" % (evt_tree.emtfUnpTrack_mode[j], evt_tree.emtfUnpTrack_eta[j], evt_tree.emtfUnpTrack_phi[j], evt_tree.emtfUnpTrack_pt[j]))

print("Total Tracks: " + str(num_tracks) + ", Mismatched Tracks: " + str(mismatched_tracks))
