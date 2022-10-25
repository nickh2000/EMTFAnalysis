import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *

'''This script is generally used to combine multiple generated plots into a single PDF'''
def cms_latex():
  cms_label = TLatex()
  cms_label.SetTextSize(0.04)
  cms_label.DrawLatexNDC(0.1, 0.92, "#bf{ #font[22]{CMS} #font[72]{Efficiency Studies}}");
  return cms_label


def head():
  header = TLatex()
  header.SetTextSize(0.03)
  header.DrawLatexNDC(0.53, 0.92, "#sqrt{s} = 13.6 TeV, Run 3 Data");

  
  return header


#plots indexed by pt-range, eta-region, and mode
#pt_ranges = ["10_20", "0_5", "5_10", "20_50", "0_50", "anypt", "0_20"]
pt_ranges = ["20_50"]
eta_ranges = ["1_1", "1_2", "2_2"]
#modes = [range(16), [11, 13, 14, 15], [7, 9, 10, 11, 13, 14, 15], [7], [9], [10], [11], [13], [14], [15]]
modes = [[15]]

infile_ideal = TFile("plots/alignment/alignment_study.root")
infile_R2 = TFile("plots/alignment/alignment_study_unpacked.root")

#index by pt range, eta-region, endcap, and mode
for pt in pt_ranges:
  for eta in eta_ranges:
    for endcap in range(2):
      for mode in modes:

        #Access plots by plots name
        PLOT_NAME = "pt_resolution_histo_e" + str(endcap) + "_" + pt + "_" + eta + "_" + str(mode) + "_"


        #Get plots from Ideal and R2Geometry files
        reso_histo_R2 = infile_R2.Get(PLOT_NAME + ";1")

        
        reso_histo_ideal = infile_ideal.Get(PLOT_NAME + ";1")
        

        R2_integral = reso_histo_R2.Integral()
        ideal_integral = reso_histo_ideal.Integral()

        if R2_integral > 0:
          reso_histo_R2.Scale(1/R2_integral)
        if ideal_integral > 0:
          reso_histo_ideal.Scale(1/ideal_integral)

        #Create PDF canvas for overlaying ideal and Run2 pt-reso plots
        canvas = TCanvas(reso_histo_R2.GetName() , reso_histo_R2.GetName(), 700,700)

        leg = TLegend(0.63, 0.7, 0.83, 0.78)
        leg.AddEntry(reso_histo_R2, 'Run 2 Geometry', 'p')

        leg.AddEntry(reso_histo_ideal, 'Run 3 Geometry', 'p')
        leg.SetMargin(0.5)
        leg.SetFillStyle(0)

        gStyle.SetOptStat(0)
        gStyle.SetLegendBorderSize(0)
        gStyle.SetLegendTextSize(0.018);

        gPad.SetGridy(1)
        gPad.SetGridx(1)

        #Configure draw style of these plots
        reso_histo_R2.SetLineWidth(0)
        reso_histo_ideal.SetLineWidth(0)

        reso_histo_R2.SetMarkerColor(kRed)
        reso_histo_ideal.SetMarkerColor(4)

        reso_histo_R2.SetMarkerStyle(20)
        reso_histo_ideal.SetMarkerStyle(21)

        reso_histo_R2.SetMarkerSize(1.3)

        #Set height of plot depending on maximum data-point
        maximum = max(reso_histo_ideal.GetMaximum(), reso_histo_R2.GetMaximum())
        reso_histo_R2.GetYaxis().SetRangeUser(0,maximum * 1.2)
        reso_histo_ideal.GetYaxis().SetRangeUser(0, maximum * 1.2)

        


        

        #Draw plots to canvas
        reso_histo_R2.Draw("g P same")
        reso_histo_ideal.Draw("g P same")
        

        #X and Y axis labels
        leg.Draw("same")
        reso_histo_R2.GetXaxis().SetTitle("(P_{T}^{L1} * q^{L1} * q^{reco}  - P_{T}^{reco}) / P_{T}^{reco}");
        reso_histo_R2.GetYaxis().SetTitle("# Muons");
        reso_histo_R2.GetYaxis().SetLabelSize(.03)

        #Set height of plot depending on maximum data-point
        maximum = max(reso_histo_ideal.GetMaximum(), reso_histo_R2.GetMaximum())
        reso_histo_R2.GetYaxis().SetRangeUser(0,maximum * 1.2)
        reso_histo_ideal.GetYaxis().SetRangeUser(0, maximum * 1.2)


        #Get gaussian fits and means of plots
        # gauss_fit_R2 = reso_histo_R2.Fit('gaus', "S Q N", reso_histo_R2.GetName() + "_gauss", -.3, .3)
        # if int(gauss_fit_R2) == -1: R2_mean = 0
        # else: R2_mean = gauss_fit_R2.Parameters()[1]


        # gauss_fit_ideal = reso_histo_ideal.Fit('gaus', "N S Q", reso_histo_ideal.GetName() + "_gauss", -.3, .3)      
        # if int(gauss_fit_ideal) == -1: ideal_mean = 0
        # else:
        #   ideal_mean = gauss_fit_ideal.Parameters()[1]



        cms_label = cms_latex()
        header = head()


        #Label plots based on pt range and eta range
        if eta == "1_1":
          if endcap == 0:
            lower = -1.6
            upper = -1.2
          elif endcap == 1:
            lower = 1.2
            upper = 1.6
        elif eta == "1_2":
          if endcap == 0:
            lower = -2.1
            upper = -1.6
          elif endcap == 1:
            lower = 1.6
            upper = 2.1
        elif eta == "2_2":
          if endcap == 0:
            lower = -2.4
            upper = -2.1
          elif endcap == 1:
            lower = 2.1
            upper = 2.4

        header.DrawLatexNDC(0.4, 0.83, "%s Gev < P_{T} < %s GeV, %s < #eta < %s" % (pt[:2], pt[-2:], lower, upper));

        # stat_box = stats()

        # stat_box.DrawLatexNDC(0.2, 0.63, "#mu_{R2} = %6.4f" % (R2_mean));
        # stat_box.DrawLatexNDC(0.2, 0.59, "#mu_{IDEAL} = %6.4f" % (ideal_mean));

        #Save canvas as PDF

        canvas.SaveAs("plots/alignment/pdfs/%s.pdf" % (PLOT_NAME))
        del canvas