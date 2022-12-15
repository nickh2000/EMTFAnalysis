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


''' This script draw high_eta_study.root histograms to PDFs
    Overlays efficiency w.r.t. phi plots in high-eta-region for different geometries'''

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


#infile_R2 = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/high_eta_study_reverse_v3.root")
#infile_R2 = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/high_eta_study_phv3.root")
infile_R2 = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/high_eta_study_custom_vTen.root")
infile_custom = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/high_eta_study_custom_chamber.root")
#infile_custom = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/high_eta_study_custom_vFourteen.root")
nFiles = 0



pt_plot = infile_R2.Get('h_phi_err_EMTF_cut_3_1;1')
canvas = TCanvas(pt_plot.GetName() , pt_plot.GetName(), 700,700)

gPad.SetGridy(1)
gPad.SetGridx(1)

#Configure draw style of these plots
pt_plot.SetLineWidth(2)

pt_plot.SetMarkerColor(4)
pt_plot.SetLineColor(4)

pt_plot.SetMarkerStyle(21)

pt_plot.SetMarkerSize(.5)


#Draw plots to canvas
pt_plot.Draw("AP same L")

#X and Y axis labels
pt_plot.GetXaxis().SetTitle("#phi");
pt_plot.GetYaxis().SetTitle("Efficiency");
pt_plot.GetYaxis().SetLabelSize(.03)


#Title and header
cms_label = cms_latex()
header = head()

#Set labels depending on endcap
header.DrawLatexNDC(0.15, 0.53, "P_{T} > 3 GeV, 2.2 < #eta < 2.4");

#Save canvas as PDF

canvas.SaveAs("plots/high_eta_region/pdfs/%s.pdf" % ("cut_3_plot"))
del canvas

#Overlay efficiencies in high eta region w.r.t. phi
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
  leg.AddEntry(histo_R2, 'Run 3 Geometry', 'p')

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


  #Title and header
  cms_label = cms_latex()
  header = head()

  #Set labels depending on endcap
  if endcap == 1:
    header.DrawLatexNDC(0.45, 0.53, "P_{T} > 22 GeV, 2.2 < #eta < 2.4");
  else:
    header.DrawLatexNDC(0.55, 0.53, "P_{T} > 22 GeV, -2.4 < #eta < -2.2");

  #Draw overall average efficiency of geometry in high-eta region
  header.DrawLatexNDC(0.2, 0.53, "#mu_{R3} = %6.4f" % (average_R2));
  header.DrawLatexNDC(0.2, 0.49, "#mu_{Custom} = %6.4f" % (average_custom));

  #Save canvas as PDF

  canvas.SaveAs("plots/high_eta_region/pdfs/%s.pdf" % (PLOT_NAME.replace(" ", "")))
  del canvas


#Draw all other PDFs
for key in infile.GetListOfKeys():
    name = key.GetName()
    plot = infile.Get(name)

    #Placeholder labels, can change with conditionals depending on plot
    plot.GetXaxis().SetTitle("Phi");
    plot.GetYaxis().SetTitle("L1T Efficiency");
    plot.GetYaxis().SetRangeUser(0.00001,1.2);
    plot.GetXaxis().SetRangeUser(-3.14159,3.14159);

    canvas = TCanvas(name, name, 700,700)


    #different draw settings for error plots
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

    
    #Write PDF's for error efficiency plots
    canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/high_eta_region/pdfs/" + name + ".pdf")


del infile