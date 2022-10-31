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


''' Script combines tempory batch output files into one final product file, for parallelization
    Adds files together, and also creates efficiency quotient plots by dividing num and den plots
    '''

#Take in arguments
#--output: file to save results too, should have the same name as the temporary input files (but exclude the index number in the temporary file's name)
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile", required = True)
parser.add_argument("-n", "--num_jobs", required = False)
args = parser.parse_args()

## Configuration settings
MAX_FILE = -1    ## Maximum number of input files (use "-1" for unlimited)
MAX_EVT  = -1       ## Maximum number of events to process
PRT_EVT  = 10000     ## Print every Nth event


#Directory to get tmp files from 
folder ="/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/tmp/"
nFiles = 0


#use num jobs argument if exists
if args.num_jobs:
    NUM_JOBS = int(args.num_jobs)
else: NUM_JOBS = 50

histos = {}
for file in Popen(['ls', folder], stdout=PIPE).communicate()[0].split():
    #only take files in tmp that contain the output file's name
    for i in range(NUM_JOBS):
        if args.outfile.split("/")[-1][:-5] + str(i) + ".root" in file: break
    else: continue

    file_name = "%s%s" % (folder, file)
    nFiles   += 1
    print ('* Loading file #%s: %s' % (nFiles, file))

    #if file is corrupt, skip it
    try: f = TFile(file_name)
    except:
        print("Skipping file %s" % file)
        continue

    #Add histograms together and store to histos directory
    for key in f.GetListOfKeys():
        name = key.GetName()
        new_histo = f.Get(name)
        if 'TGraphAsymmErrors' in str(type(new_histo)): continue
        new_histo.SetDirectory(gROOT)
        if name not in histos:
            histos[name] = new_histo
        else:
            histos[name].Add(new_histo)

#output file to which we save our sums
outfile = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/%s" % (args.outfile), "recreate")




#but we will also divide plots containing "num" or "den" in their name
for k in sorted(list(histos.keys())):
    histos[k].Write()
    if 'num' in k:
        template_string = k.replace("num", "")

        for k2, v2 in histos.items():
            if 'den' in k2:
                #denominator plots should be same as numerator plots save the keyword "den"
                if k2.replace('den', "") != template_string: continue


                #create efficiency and error plots by dividing numerator plot by denominator plot
                eff_name = k2.replace('den', 'eff')
                err_graph = TGraphAsymmErrors(histos[k], v2)
                err_graph.SetName(eff_name.replace("eff", "err"))
                histos[k].Divide(v2)
                histos[k].SetName(eff_name)
                
                if 'TH2' in histos[k].Class_Name():
                    histos[k].Draw("colz")

                #Write efficiency plot and error plot to output file
                histos[k].Write()
                err_graph.Write()

                break

del outfile