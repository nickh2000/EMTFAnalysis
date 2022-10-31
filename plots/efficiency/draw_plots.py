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


''' This script exists to draw PDFs for analyzing EMTF efficiency
    Divides efficiency w.r.t. phi&eta plots in high-eta-region for different geometries
    '''

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


infile_R2 = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/efficiency/eff_EMTF_custom.root")
infile_R3 = TFile("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/efficiency/eff_EMTF_unpacked.root")
nFiles = 0


gStyle.SetOptStat(0)
gStyle.SetLegendBorderSize(0)
gStyle.SetLegendTextSize(0.018);



#Compare total EMTF efficiency of different geometries
eta_num1 = infile_R2.Get("h_eta_num_EMTF;1")
eta_den1 = infile_R2.Get("h_eta_den_EMTF;1")
average1 = eta_num1.Integral() / eta_den1.Integral() * 2

eta_num2 = infile_R3.Get("h_eta_num_EMTF;1")
eta_den2 = infile_R3.Get("h_eta_den_EMTF;1")
average2 = eta_num2.Integral() / eta_den2.Integral() * 2

#Find ratio between geometry efficiencies
print(average1 / average2)




#Now we will draw eta-phi plots of new geometry divided by old geometry
eta_phi_custom = infile_R2.Get("h_eta_phi_eff_EMTF;1")
eta_phi_R3 = infile_R3.Get("h_eta_phi_eff_EMTF;1")

custom_integral_1 = eta_phi_custom.Integral(0, 13, 0, 30)
custom_integral_2 = eta_phi_custom.Integral(37, 50, 0, 30)
custom_integral = custom_integral_1 + custom_integral_2

R3_integral_1 = eta_phi_R3.Integral(0, 13, 0, 29)
R3_integral_2 = eta_phi_R3.Integral(37, 50, 0, 50)
R3_integral = R3_integral_1 + R3_integral_2


canvas = TCanvas(eta_phi_custom.GetName() , eta_phi_custom.GetName(), 700,700)

eta_phi_custom.Divide(eta_phi_R3)

eta_phi_custom.Draw("colz")

cms_label = cms_latex()
header = head()

R3_average = R3_integral / (26.*30.)
custom_average = custom_integral / (26.*30.)
header.DrawLatexNDC(0.2, 0.63, "#mu_{R3} = %6.4f" % (R3_average));
header.DrawLatexNDC(0.2, 0.59, "#mu_{Custom} = %6.4f" % (custom_average));

canvas.SaveAs("plots/efficiency/pdfs/%s.pdf" % ("eta_phi_compare"))


#Access plots by plots name
PLOT_NAME = "h_pt_err_EMTF_large"


#Get plots from R3 and R2 Geometry files
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



#Draw all other plots to pdf
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
    
    #Write PDF's for error efficiency plots
    canvas.SaveAs("/afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/plots/efficiency/pdfs/" + name + ".pdf")


del infile