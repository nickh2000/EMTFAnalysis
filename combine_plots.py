import sys
import math
from subprocess import Popen,PIPE
from ROOT import *
import numpy as np
from array import *
import Helper as h




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

def stats():
  stats = TLatex()
  stats.SetTextSize(0.025)
  
  return stats


'''This section devides unpacked efficiencies by packed efficiencies'''
infile = TFile.Open( 'plots/L1T_eff.root',"READ")
packed = infile.Get('h_eta_phi_trg_EMTF_eff;')


canvas = canvas = TCanvas(packed.GetName() , packed.GetName(), 700,700)
packed.Draw("colz")
gStyle.SetOptStat(0)
canvas.SaveAs("plots/phi_eta_eff_run2.pdf")



if not packed :
    print(" Failed to get data histogram ")


packed.SetDirectory(0)
del infile

infile = TFile.Open( 'plots/L1_eff_run3.root',"READ")
unpacked = infile.Get('h_eta_phi_trg_EMTF_eff;1')

canvas = canvas = TCanvas(unpacked.GetName() , unpacked.GetName(), 700,700)
unpacked.Draw("colz")
gStyle.SetOptStat(0)
canvas.SaveAs("plots/phi_eta_eff_run3.pdf")

unpacked.SetDirectory(0)
del infile


print((unpacked.Integral(38, 50, 1, 30) + unpacked.Integral(1, 12, 1, 30)) / 24 / 30 / ((packed.Integral(38, 50, 1, 30) + packed.Integral(1, 12, 1, 30)) / 24 / 30))

out_file = TFile('plots/eff_divided.root', 'recreate')
unpacked = unpacked.Clone()
unpacked.Divide(packed)
unpacked.SetName('unpacked-packed')

canvas = TCanvas(unpacked.GetName() , unpacked.GetName(), 700,700)
unpacked.Draw("colz")

unpacked.SetLineWidth(1)
unpacked.SetLineColor(1)

cms_label =cms_latex()
header = head()

gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018)
gStyle.SetOptStat(0)

canvas.SaveAs("plots/" + unpacked.GetName() + ".pdf")

unpacked.Write()

print((unpacked.Integral(38, 50, 1, 30) + unpacked.Integral(1, 12, 1, 30)) / 24 / 30)
del out_file




#plots indexed by pt-range, eta-region, and mode
pt_ranges = ["10_20", "0_5", "5_10", "20_50", "0_50", "anypt", "0_20"]
eta_ranges = ["1_1", "1_2", "2_2"]
modes = [range(16), [11, 13, 14, 15], [7, 9, 10, 11, 13, 14, 15], [7], [9], [10], [11], [13], [14], [15]]


#Source plots if we want zerobias dataset
ZERO_BIAS = False
RUN3 = True
UNP = True

if ZERO_BIAS:
  infile_ideal = TFile.Open( 'plots/dphi_zerobias_ideal.root',"READ")
  infile_R2 = TFile.Open( 'plots/dphi_zerobias.root',"READ")
else:
  infile_ideal = TFile.Open( 'plots/dphi_ideal.root',"READ")
  infile_R2 = TFile.Open( 'plots/dphi.root',"READ")




if RUN3:
  infile_ideal = TFile.Open( 'plots/dphi_run3.root',"READ")
  #Use unpacked data instead of R2 data
  if UNP:
    infile_R2 = TFile.Open( 'plots/dphi_run3_unp.root',"READ")


#Get rid of stats box and confiture legend
gStyle.SetOptStat(0)
gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018)
      


#index by pt range, eta-region, endcap, and mode
for pt in pt_ranges:
  for eta in eta_ranges:
    for endcap in range(2):
      for mode in modes:

        #Access plots by plots name
        PLOT_NAME = "pt_resolution_histo_e" + str(endcap) + "_" + pt + "_" + eta + "_" + str(mode) + "_"


        #Get plots from Ideal and R2Geometry files
        reso_histo_R2 = infile_R2.Get(PLOT_NAME + ";1")

        if not RUN3:
          reso_histo_ideal = infile_ideal.Get(PLOT_NAME + "ideal;1")
        else:   reso_histo_ideal = infile_ideal.Get(PLOT_NAME + ";1")
        

        R2_integral = reso_histo_R2.Integral()
        ideal_integral = reso_histo_ideal.Integral()

        if R2_integral > 0:
          reso_histo_R2.Scale(1/R2_integral)
        if ideal_integral > 0:
          reso_histo_ideal.Scale(1/ideal_integral)

        #Create PDF canvas for overlaying ideal and Run2 pt-reso plots
        canvas = TCanvas(reso_histo_R2.GetName() , reso_histo_R2.GetName(), 700,700)

        leg = TLegend(0.63, 0.7, 0.83, 0.78)
        leg.AddEntry(reso_histo_R2, 'Run2 Geometry', 'p')

        if not RUN3: leg.AddEntry(reso_histo_ideal, 'Ideal Geometry', 'p')
        else: leg.AddEntry(reso_histo_ideal, 'Run 3 Geometry', 'p')
        leg.SetMargin(0.5)
        leg.SetFillStyle(0)

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
        reso_histo_R2.Draw("P same")
        reso_histo_ideal.Draw("P same")
        

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

        if ZERO_BIAS:
          canvas.SaveAs("plots/resolutions/zerobias_%s.pdf" % (PLOT_NAME))
        elif RUN3:
          canvas.SaveAs("plots/resolutions/RUN3vsRUN2/%s.pdf" % (PLOT_NAME))
        del canvas

  