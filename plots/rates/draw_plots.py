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


infile = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/rates/rate_calculation_unpacked.root")
nFiles = 0

nums = np.zeros(3)
denom = 0

#Create PDF canvas for overlaying ideal and Run2 pt-distribution plots
canvas = TCanvas('pt_1D' , 'pt_1D', 700,700)

gStyle.SetOptStat(0)
gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018);

gPad.SetGridy(1)
gPad.SetGridx(1)
gPad.SetLogy(1)

leg = TLegend(0.63, 0.7, 0.83, 0.78)
leg.SetMargin(0.5)
leg.SetFillStyle(0)
for offset in range(-2, 1):

  #1D pt distribution plot
  plot = infile.Get('pt_1D_%d' % offset)
  denominator = infile.Get('bx_den')

  #scale by the number of tracks in that bx
  plot.Scale( 1 / denominator.GetBinContent(offset + 3))

  leg.AddEntry(plot, 'BX %d' % offset, 'p')

  
  #Configure draw style of these plots
  plot.SetLineWidth(0)

  if offset == -2: plot.SetMarkerColor(kRed)
  if offset == -1: plot.SetMarkerColor(kBlue)
  if offset == 0: plot.SetMarkerColor(kGreen)

  plot.SetMarkerStyle(20)

  plot.SetMarkerSize(1.3)

  #Draw plots to canvas

  plot.Draw("g P same")
  
  #X and Y axis labels
  
  plot.GetXaxis().SetTitle("PT");
  plot.GetYaxis().SetTitle("# Tracks");
  plot.GetYaxis().SetLabelSize(.03)

cms_label = cms_latex()
header = head()

leg.Draw("same")
canvas.SaveAs("plots/rates/pdfs/pt_1D.pdf" % ())
del canvas



#Create PDF canvas for overlaying ideal and Run2 mode-distribution plots
canvas = TCanvas('mode_1D' , 'mode_1D', 700,700)

gStyle.SetOptStat(0)
gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018);

gPad.SetGridy(1)
gPad.SetGridx(1)

leg = TLegend(0.13, 0.7, 0.33, 0.78)
leg.SetMargin(0.5)
leg.SetFillStyle(0)
for offset in range(-2, 1):
  plot = infile.Get('mode_1D_%d' % offset)
  denominator = infile.Get('bx_den')

  #normalize plot by number of tracks in that BX
  plot.Scale(1 / denominator.GetBinContent(offset + 3))

  leg.AddEntry(plot, 'BX %d' % offset, 'p')

  
  #Configure draw style of these plots
  plot.SetLineWidth(0)

  if offset == -2: plot.SetMarkerColor(kRed)
  if offset == -1: plot.SetMarkerColor(kBlue)
  if offset == 0: plot.SetMarkerColor(kGreen)

  plot.SetMarkerStyle(20)

  plot.SetMarkerSize(1.3)

  #Draw plots to canvas

  plot.Draw("g P same")
  
  #X and Y axis labels
  
  plot.GetXaxis().SetTitle("Mode");
  plot.GetYaxis().SetTitle("# Tracks");
  plot.GetYaxis().SetLabelSize(.03)

cms_label = cms_latex()
header = head()

leg.Draw("same")
canvas.SaveAs("plots/rates/pdfs/mode_1D.pdf" % ())
del canvas


#Draw all the other plots to a PDF
for key in infile.GetListOfKeys():
    name = key.GetName()
    plot = infile.Get(name)

    if 'pt_1D' in name: continue

    #get incidence of all tracks across all BXs
    if 'den' in name: denom = plot.Integral()

    #Get incidence of tracks per BX
    if 'num' in name:
      offset = int(name.split(';')[0].split('_')[-1])
      nums[offset + 2] = plot.Integral()

    #Normalize plots
    if 'high_eta_bx' in name:
      denominator = infile.Get('bx_den')
      plot.Divide(denominator)
    if 'bx_rpc' in name: 
      denominator = infile.Get('bx_den')
      plot.Divide(denominator)
    if 'high_eta_mode' in name:
      offset = int(name.split(';')[0].split('_')[-1])
      denominator = infile.Get("mode_den_%d" % offset)
      plot.Divide(denominator)

    #Draw these to PDF
    canvas = TCanvas(plot.GetName() , plot.GetName(), 700,700)

    gStyle.SetOptStat(0)

    plot.Draw("colz")
    plot.GetXaxis().SetTitle("#eta")
    plot.GetYaxis().SetTitle("#phi")
    cms_label = cms_latex()
    header = head()

    canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/rates/pdfs/" + name + ".pdf")

#Print rates of different BXs
print("BX = 0 Rate: " + str(nums[2] / denom))
print("BX = -1 Rate: " + str(nums[1] / denom))
print("BX = -2 Rate: " + str(nums[0] / denom))

del infile