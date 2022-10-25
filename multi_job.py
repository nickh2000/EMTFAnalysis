from cmath import tan
from re import L
import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
import Helper as h
import os
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
parser.add_argument("-o", "--outfile", required = True)
parser.add_argument("-n", "--num_jobs", required = False)
args = parser.parse_args()

## Configuration settings
MAX_FILE = -1    ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event


#Specify which datasets to use
RUN3 = True
ZERO_BIAS = False


folder ="/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/"
nFiles = 0

if args.num_jobs:
    NUM_JOBS = int(args.num_jobs)
else: NUM_JOBS = 50

histos = {}
for file in Popen(['ls', folder], stdout=PIPE).communicate()[0].split():
    for i in range(NUM_JOBS):
        if args.outfile.split("/")[-1][:-5] + str(i) + ".root" in file: break
    else: continue

    file_name = "%s%s" % (folder, file)
    nFiles   += 1
    print ('* Loading file #%s: %s' % (nFiles, file))
    try: f = TFile(file_name)
    except:
        print("Skipping file %s" % file)
        continue
    for key in f.GetListOfKeys():
        name = key.GetName()
        new_histo = f.Get(name)
        if 'TGraphAsymmErrors' in str(type(new_histo)): continue
        new_histo.SetDirectory(gROOT)
        if name not in histos:
            histos[name] = new_histo
        else:
            histos[name].Add(new_histo)

outfile = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/%s" % (args.outfile), "recreate")





for k in sorted(list(histos.keys())):
    print(k)
    histos[k].Write()
    if 'num' in k:
        template_string = k.replace("num", "")

        for k2, v2 in histos.items():
            if 'den' in k2:
                if k2.replace('den', "") != template_string: continue

                eff_name = k2.replace('den', 'eff')
                eff_graph = TGraphAsymmErrors(histos[k], v2)
                eff_graph.SetName(eff_name.replace("eff", "err"))
                histos[k].Divide(v2)
                print(k)
                histos[k].SetName(eff_name)
                

                print(histos[k].Class_Name())
                if 'TH2' in histos[k].Class_Name():
                    histos[k].Draw("colz")
                histos[k].Write()
                eff_graph.Write()

                break



del outfile

# deleted_folder = "/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/"

# for filename in os.listdir(deleted_folder):
#     file_path = os.path.join(deleted_folder, filename)
#     try:
#         if os.path.isfile(file_path) or os.path.islink(file_path):
#             os.unlink(file_path)
#         elif os.path.isdir(file_path):
#             shutil.rmtree(file_path)
#     except Exception as e:
#         print('Failed to delete %s. Reason: %s' % (file_path, e))
    

# if __name__ == "__main__":
#     folder ="/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/"

#     files = []

#     histos = {}
#     for file in Popen(['ls', folder], stdout=PIPE).communicate()[0].split():
#         for i in range(50):
#             if args.infile[:-3] + str(i) + ".root" in file: break
#         else: continue
#         files.append(file)

#         file_name = "%s%s" % (folder, file)
#         nFiles   += 1
#         print ('* Loading file #%s: %s' % (nFiles, file))

#         del f

#     Hadd('%s' % (args.outfile[:-3] + ".root", files))

