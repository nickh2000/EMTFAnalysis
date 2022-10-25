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


infile_R2 = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/high_eta_study_phv3.root")
infile_custom = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/high_eta_study_ref_2.root")
nFiles = 0


for endcap in range(2):

  gStyle.SetOptStat(0)
  gStyle.SetLegendBorderSize(0)
  gStyle.SetLegendTextSize(0.018)


  #Access plots by plots name
  PLOT_NAME_NUM = "h_phi_num_EMTF_mode_[11, 13, 14, 15]_%d" % endcap
  PLOT_NAME_DEN = "h_phi_den_EMTF_mode_[11, 13, 14, 15]_%d" % endcap
  PLOT_NAME = "h_phi_err_EMTF_mode_[11, 13, 14, 15]_%d" % endcap


  #Get plots from Ideal and R2Geometry files
  histo_R2 = infile_R2.Get(PLOT_NAME + ";1")
  histo_R2_NUM = infile_R2.Get(PLOT_NAME_NUM + ";1")
  histo_R2_DEN = infile_R2.Get(PLOT_NAME_DEN + ";1")
  average_R2 = histo_R2_NUM.Integral() / histo_R2_DEN.Integral()

  histo_custom = infile_custom.Get(PLOT_NAME + ";1")
  histo_custom_NUM = infile_custom.Get(PLOT_NAME_NUM + ";1")
  histo_custom_DEN = infile_custom.Get(PLOT_NAME_DEN + ";1")
  average_custom = histo_custom_NUM.Integral() / histo_custom_DEN.Integral()


  #Create PDF canvas for overlaying ideal and Run2 pt-reso plots
  canvas = TCanvas(histo_R2.GetName() , histo_R2.GetName(), 700,700)

  leg = TLegend(0.63, 0.3, 0.83, 0.38)
  leg.AddEntry(histo_R2, 'Run 2 Geometry', 'p')

  leg.AddEntry(histo_custom, 'Custom Geometry', 'p')
  leg.SetMargin(0.5)
  leg.SetFillStyle(0)

  

  gPad.SetGridy(1)
  gPad.SetGridx(1)

  #Configure draw style of these plots
  histo_R2.SetLineWidth(2)
  histo_custom.SetLineWidth(2)

  histo_R2.SetMarkerColor(kRed)
  histo_R2.SetLineColor(kRed)
  histo_custom.SetMarkerColor(4)
  histo_custom.SetLineColor(4)

  histo_R2.SetMarkerStyle(20)
  histo_custom.SetMarkerStyle(21)

  histo_R2.SetMarkerSize(.5)
  histo_custom.SetMarkerSize(.5)


  #Draw plots to canvas
  histo_R2.Draw("AP same L")
  histo_custom.Draw("p l same")


  #X and Y axis labels
  leg.Draw("same")
  histo_R2.GetXaxis().SetTitle("#phi");
  histo_R2.GetYaxis().SetTitle("Efficiency");
  histo_R2.GetYaxis().SetLabelSize(.03)

  #Set height of plot depending on maximum data-point
  maximum = max(histo_custom.GetMaximum(), histo_R2.GetMaximum())
  histo_R2.GetYaxis().SetRangeUser(0,maximum * 1.2)
  histo_custom.GetYaxis().SetRangeUser(0, maximum * 1.2)


  # histo_R2.GetXaxis().SetRangeUser(-3.14159, 3.14159)
  # histo_custom.GetXaxis().SetRangeUser(-3.14159, 3.14159)

  cms_label = cms_latex()
  header = head()


  header.DrawLatexNDC(0.4, 0.83, "P_{T} < 22 GeV, 2.2 < #eta < 2.4");


  header.DrawLatexNDC(0.2, 0.63, "#mu_{R2} = %6.4f" % (average_R2));
  header.DrawLatexNDC(0.2, 0.59, "#mu_{Custom} = %6.4f" % (average_custom));

  #Save canvas as PDF

  canvas.SaveAs("plots/high_eta_region/pdfs/%s.pdf" % (PLOT_NAME.replace(" ", "")))
  del canvas

for key in infile.GetListOfKeys():
    name = key.GetName()
    plot = infile.Get(name)

    plot.GetXaxis().SetTitle("Phi");
    plot.GetYaxis().SetTitle("L1T Efficiency");
    plot.GetYaxis().SetRangeUser(0.00001,1.2);
    plot.GetXaxis().SetRangeUser(-3.14159,3.14159);

    canvas = TCanvas(name, name, 700,700)

    if 'err' in name: plot.Draw("A P C")
    else: plot.Draw("g P C")

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
    canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/pdfs/" + name + ".pdf")


del infile