from cmath import tan
from re import L
import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
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


infile_R3 = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/efficiency/eff_EMTF.root")
infile_R2 = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/efficiency/eff_EMTF_unpacked.root")
nFiles = 0

#Access plots by plots name
PLOT_NAME = "h_pt_err_EMTF_large"


#Get plots from Ideal and R2Geometry files
reso_histo_R2 = infile_R2.Get(PLOT_NAME + ";1")


reso_histo_R3 = infile_R3.Get(PLOT_NAME + ";1")

#Create PDF canvas for overlaying ideal and Run2 pt-reso plots
canvas = TCanvas(reso_histo_R2.GetName() , reso_histo_R2.GetName(), 700,700)

leg = TLegend(0.63, 0.8, 0.83, 0.88)
leg.AddEntry(reso_histo_R2, 'Run 2 Geometry', 'p')

leg.AddEntry(reso_histo_R3, 'Run 3 Geometry', 'p')
leg.SetMargin(0.5)
leg.SetFillStyle(0)
leg.SetBorderSize(0)

gStyle.SetOptStat(0)
gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018);

gPad.SetGridy(1)
gPad.SetGridx(1)

#Configure draw style of these plots
reso_histo_R2.SetLineWidth(2)
reso_histo_R3.SetLineWidth(2)

reso_histo_R2.SetLineColor(kRed)
reso_histo_R3.SetLineColor(4)

reso_histo_R2.SetMarkerStyle(20)
reso_histo_R3.SetMarkerStyle(21)

reso_histo_R2.SetMarkerColor(kRed)
reso_histo_R3.SetMarkerColor(4)

reso_histo_R2.SetMarkerSize(.8)
reso_histo_R3.SetMarkerSize(.8)

#Draw plots to canvas
reso_histo_R2.Draw("APL")
reso_histo_R3.Draw("p same")


#X and Y axis labels
leg.Draw("same")
reso_histo_R2.GetXaxis().SetTitle("P_{T}");
reso_histo_R2.GetXaxis().SetRangeUser(0, 1000);
reso_histo_R2.GetYaxis().SetTitle("Efficiency");
reso_histo_R2.GetYaxis().SetLabelSize(.03)

reso_histo_R2.GetYaxis().SetRangeUser(0, 1.2)
reso_histo_R3.GetYaxis().SetRangeUser(0, 1.2)


cms_label = cms_latex()
header = head()

#Save canvas as PDF

canvas.SaveAs("plots/efficiency/pdfs/%s.pdf" % (PLOT_NAME))
del canvas
del infile_R3
del infile_R2


infile = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/efficiency/eff_EMTF.root")

for key in infile.GetListOfKeys():
    name = key.GetName()
    plot = infile.Get(name)

    if 'err' in name: continue
    else: plot.Draw("g P C")
    plot.GetXaxis().SetTitle("Phi");
    plot.GetYaxis().SetTitle("L1T Efficiency");
    plot.GetYaxis().SetRangeUser(0.00001,1.2);
    plot.GetXaxis().SetRangeUser(-3.14159,3.14159);

    canvas = TCanvas(name, name, 700,700)

    

    plot.SetLineWidth(1)
    plot.SetLineColor(1)

    cms_label =cms_latex()
    header = head()

    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018);
    gStyle.SetOptStat(0)

    leg =TLegend(0.4,0.8,0.88,0.88);

    # leg.AddEntry(plot, "L1 pT (prompt) > 22 GeV, eta < -2.2")
    # leg.Draw("same P")
    
    #Write PDF's for error efficiency plots
    canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/efficiency/pdfs/" + name + ".pdf")


del infile