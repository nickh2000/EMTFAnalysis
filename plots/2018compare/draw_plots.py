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


infile = TFile("plots/2018compare/2018compare.root")


def compare(property):
    nFiles = 0

    #Overlay efficiencies in high eta region w.r.t. phi

    gStyle.SetOptStat(0)
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)

    #Get plots from Ideal and R2Geometry files

    PLOT_NAME = "%s_EMTF_2022" % property
    PLOT_NAME2 = "%s_EMTF_2018" % property

    histo= infile.Get(PLOT_NAME)
    histo.Scale(1./histo.Integral())
    histo2 = infile.Get(PLOT_NAME2)
    histo2.Scale(1./histo2.Integral())

    #Create PDF canvas for overlaying ideal and Run2 pt-reso plots
    canvas = TCanvas(histo.GetName() , histo.GetName(), 700,700)



    leg = TLegend(0.63, 0.5, 0.83, 0.78)
    leg.AddEntry(histo, '2022', 'p')
    leg.AddEntry(histo2, '2018', 'p')
    leg.SetMargin(0.5)
    leg.SetFillStyle(0)
    
    gPad.SetGridy(1)
    gPad.SetGridx(1)

    

    histo2.SetMarkerColor(kBlue)
    histo2.SetLineWidth(2)
    histo2.SetLineColor(kBlue)
    histo2.SetMarkerStyle(21)
    histo2.SetMarkerSize(.5)

    #Configure draw style of these plots
    histo.SetMarkerColor(kRed)
    histo.SetLineColor(kRed)
    histo.SetMarkerStyle(20)
    histo.SetLineWidth(2)
    histo.SetMarkerSize(.5)

    #Save canvas as PDF

    maximum = max(histo.GetMaximum(), histo2.GetMaximum())

    histo.GetYaxis().SetRangeUser(10e-7, 1.1 * maximum)
    histo.GetXaxis().SetTitle(property)


    if 'pt' in property:
      gPad.SetLogy()


    histo.Draw("hist P L")
    histo2.Draw("hist P L same")
    leg.Draw("same")


    canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/2018compare/pdfs/compare_%s.pdf" % property)
    del canvas

def divide_plots(property):
    #Overlay efficiencies in high eta region w.r.t. phi

    gStyle.SetOptStat(0)
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)

    #Get plots from Ideal and R2Geometry files

    PLOT_NAME = "%s_EMTF_2022" % property
    PLOT_NAME2 = "%s_EMTF_2018" % property

    histo= infile.Get(PLOT_NAME)
    #histo.Scale(1./histo.Integral())
    histo2 = infile.Get(PLOT_NAME2)
    #histo2.Scale(1./histo2.Integral())


    histo.Divide(histo2)
    histo.Scale(1./histo.Integral())


    #Create PDF canvas for overlaying ideal and Run2 pt-reso plots
    canvas = TCanvas(histo.GetName() , histo.GetName(), 700,700)

    #histo.GetYaxis().SetRangeUser(0, 60)
    histo.SetTitle(property)

    if 'pt' in property:
      histo.GetYaxis().SetRangeUser(0, 60)
    gPad.SetLogz()

    histo.Draw("colz")
    outfile = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/2018compare/more.root", "recreate")
    histo.Write()
    del outfile

    canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/2018compare/pdfs/compare_%s.pdf" % property)
    del canvas

compare('pt')
compare('pt_2')
#qual >= 12
compare('pt_2_0')
#qual >= 8
compare('pt_2_1')
#qual >= 4
compare('pt_2_2')

#1.2 < eta < 1.6, qual >= 12
compare('pt_2_0_0')
#1.2 < eta < 1.6, qual >= 8
compare('pt_2_0_1')
#1.2 < eta < 1.6, qual >= 4
compare('pt_2_0_2')

#1.6 < eta < 2.1, qual >= 12
compare('pt_2_1_0')
#1.6 < eta < 2.1, qual >= 8
compare('pt_2_1_1')
#1.6 < eta < 2.1, qual >= 4
compare('pt_2_1_2')

#2.1 < eta < 2.5, qual >= 12
compare('pt_2_2_0')
#2.1 < eta < 2.5, qual >= 8
compare('pt_2_2_1')
#2.1 < eta < 2.5, qual >= 4
compare('pt_2_2_2')

compare('eta')
compare('phi')
compare('rate-phi')
compare('rate-eta')
divide_plots('pt-phi')
divide_plots('pt-eta')
divide_plots('rate-eta-phi')


# canvas = TCanvas(histo.GetName() , histo.GetName(), 700,700)
# histo.SetTitle("Run 2022 / Run 2018")
# gPad.SetGridy(1)
# gPad.SetGridx(1)
# gPad.SetLogy()
# histo.Divide(histo2)
# histo.Draw("hist C")
# canvas.SaveAs("pdfs/dividePT.pdf")

del infile