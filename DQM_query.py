import os
from requests_futures.sessions import FuturesSession
import requests
import io
import cgi
from ROOT import *
from bs4 import BeautifulSoup
import json
import runregistry

#Define Headers and Label###############################################
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



##########Define Main Builk of Code#######################
def main():

    ##Configure session with appropriate SSL certificates
    sess = FuturesSession()
    sess.verify = "/.globus/CERN_Root_CA.crt"
    sess.cache = "/.globus"
    sess.cert = ("/.globus/usercert.pem", "/.globus/userkey.pem")

    TIMEOUT = 5

    def _fetch_dqm_rows(url, timeout=TIMEOUT):
            """Return a future of DQMRows of a DQM page at url.
            Access the array of DQMRows at _resolve(self._fetch_dqm_rows(...)).data"""

            return sess.get(url, timeout=timeout, verify = "/home/nick/.globus/CERN_Root_CA.crt", stream=True)
    


    ##Define Runs of Interest######################################
    min_run = 355100
    #min_run = 356169 # July 25
    # max_run = 357329 # August 11

    min_run = 357330 #August 12
    #max_run = 357910 #August 23, LHC shutfown
    max_run = 359763

    ########### Plot uGMT track eta's with |eta| > 2.5
    over_eta = TH1D('over_eta',  '', 50, min_run, max_run)
    over_eta_pos = TH1D('over_eta_pos',  '', 50, min_run, max_run)
    over_eta_neg = TH1D('over_eta_neg',  '', 50, min_run, max_run)


    #Plots EMTF Tracks with mode = 0
    mode_zero = TH1D('mode_zero',  '', 50, min_run, max_run)

    #Plot hardware tracks with eta bin == +/- 230
    over_etahw = TH1D('over_etahw',  '', 50, min_run, max_run)


    #Plots over_eta but indexes by all types, cosmics, and collisions respectively
    run_type_plot = [TH2D('run_type_all', '', 100, min_run, max_run, 2, 0, 2), 
                    TH1D('run_type_cosmics', '', 100, min_run, max_run), 
                    TH1D('run_type_collisions', '', 100, min_run, max_run)]

    #Plots integrals of each collisions run
    over_eta_int = TH1D('over_eta_int', '', 100, min_run, max_run)

    mismatch_plot = TH1D('mismatch_plot', '', 100, min_run, max_run)
    gem_plot = TH1D('gem_plot', '', 100, min_run, max_run)

    out_file = TFile('plots/high_eta_EMTF.root', 'recreate')

    #Get list of runs in range of type collision
    request = runregistry.get_runs(filter={
    'class': {'or': ['Collisions22']},
    'run_number':{
      'and':[
        {'>=': min_run},
        {'<=': max_run},
        ]
      }
    })

    collision_runs = {run['run_number']: run['oms_attributes']['end_time'][5:10] for run in request}
    
    #Get list of runs in range of type cosmic
    request = runregistry.get_runs(filter={
    'class': {'or': ['Commissioing22', 'Cosmics22']},
    'run_number':{
      'and':[
        {'>=': min_run},
        {'<=': max_run},
      ]
    }
  })

    cosmics_runs = {run['run_number']: run['oms_attributes']['end_time'][5:10] for run in request}




    #Go through all runs in range
    for run in range(min_run, max_run):


        #Get hisogram json for uGMT eta's
        eta_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1T/L1TStage2uGMT/ugmtMuonEta?formatted=true' % run
        

        response = _fetch_dqm_rows(eta_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})
        
        data = json.loads(elem.text)


        #Print to show progress
        if not run % 100:
            print('Processing run %d' % run)

        #move on if run is invalid
        if data['hist'] == 'unsupported type': continue

        hist = data['hist']['bins']['content']
        
        ##Print according to out-of-range eta's are positive or negative
        instances_pos = hist[-1] + data['hist']['stats']['overflow']

        instances_neg = hist[0] + data['hist']['stats']['underflow']

        over_eta.Fill(run, instances_pos + instances_neg)
        over_eta_pos.Fill(run, instances_pos)
        over_eta_neg.Fill(run, instances_neg)


        #Seperate runs by cosmics vs. collisions
        if run in cosmics_runs.keys(): 
            run_type_plot[0].Fill(run, 0, instances_neg + instances_pos)
            run_type_plot[1].Fill(run, instances_neg + instances_pos)

        elif run in collision_runs.keys(): 
          run_type_plot[0].Fill(run, 1, instances_neg + instances_pos)
          run_type_plot[2].Fill(run, instances_neg + instances_pos)

          if instances_neg + instances_pos > 0: print(run)
          over_eta_int.Fill(run, data['hist']['bins']['integral'])
        



        #Reference json with EMTF track mode hisotra
        mode_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1T/L1TStage2EMTF/emtfTrackMode?formatted=true' % run

        response = _fetch_dqm_rows(mode_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})

        data = json.loads(elem.text)

        if data['hist'] == 'unsupported type': continue

        hist = data['hist']['bins']['content']

        if run in collision_runs.keys(): 
          mode_zero.Fill(run, hist[0])

        etahw_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1T/L1TStage2EMTF/MuonCand/emtfMuonhwEta?formatted=true' % run

        response = _fetch_dqm_rows(etahw_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})

        data = json.loads(elem.text)
        
        if data['hist'] == 'unsupported type': continue
        
        hist = data['hist']['bins']['content']

        
        if run in collision_runs.keys(): 
          over_etahw.Fill(run, hist[0] + hist[-1] + data['hist']['stats']['overflow'] + data['hist']['stats']['underflow'] )
        

        mismatch_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1TEMU/L1TdeStage2EMTF/mismatchRatio?formatted=true' % run

        response = _fetch_dqm_rows(mismatch_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})

        data = json.loads(elem.text)
        
        
        if data['hist'] == 'unsupported type': continue

        hist = data['hist']['bins']['content']

        if run in collision_runs.keys(): 
          print('mismatch' + str(hist[2]))
          mismatch_plot.Fill(run, hist[2] * 100)


        
        gem_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1T/L1TStage2EMTF/gemHitOccupancy?formatted=true' % run

        response = _fetch_dqm_rows(gem_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})

        data = json.loads(elem.text)
        
        if data['hist'] == 'unsupported type': continue

        if run in collision_runs.keys(): 
          print('gem: ' + str(data['hist']['stats']['entries']))
          gem_plot.Fill(run, data['hist']['stats']['entries'])

    
    over_eta.GetXaxis().SetTitle("Run Number")
    over_eta.GetYaxis().SetTitle("Instances")
    over_eta.GetXaxis().SetLabelSize(.02)
    
    canvas = TCanvas(over_eta.GetName() , over_eta.GetName(), 900,900)
    canvas.SetLogy()
    over_eta.Draw("HIST")
    over_eta.GetXaxis().SetLabelSize(.02)
    over_eta.GetXaxis().LabelsOption('v')

    date = ""
    for run_num, new_date in collision_runs.items():
        if new_date != date:
          date = new_date
          over_eta.GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * over_eta.GetNbinsX() + 1), '%d: %s' % (run_num, date))

    over_eta.SetLineWidth(1)
    over_eta.SetLineColor(1)

    cms_label =cms_latex()
    header = head()
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)
    gStyle.SetOptStat(0)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(over_eta, "Muons with |#eta| > 2.5")
    #leg.Draw("same P")

    canvas.SaveAs("plots/over_eta.pdf")
    over_eta.Write()



    over_eta_pos.GetXaxis().SetTitle("Run Number")
    over_eta_pos.GetYaxis().SetTitle("Instances")
    over_eta_pos.GetXaxis().SetNdivisions(5)
    
    canvas = TCanvas(over_eta_pos.GetName() , over_eta_pos.GetName(), 900,900)
    canvas.SetLogy()
    over_eta_pos.Draw("HIST")

    over_eta_pos.SetLineWidth(1)
    over_eta_pos.SetLineColor(1)

    cms_label =cms_latex()
    header = head()
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(over_eta_pos, "Muons with |#eta| > 2.5")
    #leg.Draw("same P")

    canvas.SaveAs("plots/over_eta_pos.pdf")
    over_eta_pos.Write()



    over_eta_neg.GetXaxis().SetTitle("Run Number")
    over_eta_neg.GetYaxis().SetTitle("Instances")
    over_eta_neg.GetXaxis().SetNdivisions(5)
    
    canvas = TCanvas(over_eta_neg.GetName() , over_eta_neg.GetName(), 900,900)
    canvas.SetLogy()
    over_eta_neg.Draw("HIST")

    over_eta_neg.SetLineWidth(1)
    over_eta_neg.SetLineColor(1)

    cms_label =cms_latex()
    header = head()
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(over_eta_neg, "Muons with |#eta| > 2.5")
    #leg.Draw("same P")

    canvas.SaveAs("plots/over_eta_neg.pdf")
    over_eta_neg.Write()




    mode_zero.GetXaxis().SetTitle("Run Number")
    mode_zero.GetYaxis().SetTitle("Instances")
    mode_zero.GetXaxis().SetNdivisions(5)
    
    canvas = TCanvas(mode_zero.GetName() , mode_zero.GetName(), 900,900)
    canvas.SetLogy()
    mode_zero.Draw("HIST")

    date = ""
    for run_num, new_date in collision_runs.items():
        if new_date != date:
          date = new_date
          mode_zero.GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * mode_zero.GetNbinsX() + 1), '%d: %s' % (run_num, date))

    mode_zero.SetLineWidth(1)
    mode_zero.SetLineColor(1)
    mode_zero.GetXaxis().SetLabelSize(.02)
    mode_zero.GetXaxis().LabelsOption('v')

    cms_label =cms_latex()
    header = head()
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)
    gStyle.SetOptStat(0)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(mode_zero, "Muons with |#eta| > 2.5")
    #leg.Draw("same P")

    canvas.SaveAs("plots/mode_zero.pdf")
    mode_zero.Write()


    over_etahw.GetXaxis().SetTitle("Run Number")
    over_etahw.GetYaxis().SetTitle("Instances")
    canvas = TCanvas(over_etahw.GetName() , over_etahw.GetName(), 900,900)
    canvas.SetLogy()
    over_etahw.Draw("HIST")

    date = ""
    for run_num, new_date in collision_runs.items():
        if new_date != date:
          date = new_date
          over_etahw.GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * over_etahw.GetNbinsX() + 1), '%d: %s' % (run_num, date))

    over_etahw.SetLineWidth(1)
    over_etahw.SetLineColor(1)
    over_etahw.GetXaxis().SetLabelSize(.02)
    over_etahw.GetXaxis().LabelsOption('v')

    cms_label =cms_latex()
    header = head()
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)
    gStyle.SetOptStat(0)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(over_etahw, "Muons with |#eta| > 2.5")
    #leg.Draw("same P")

    canvas.SaveAs("plots/over_etahw.pdf")

    over_etahw.Write()




    gem_plot.GetXaxis().SetTitle("Run Number")
    gem_plot.GetYaxis().SetTitle("Instances")
    canvas = TCanvas(gem_plot.GetName() , gem_plot.GetName(), 900,900)
    canvas.SetLogy()
    gem_plot.Draw("HIST")

    date = ""
    for run_num, new_date in collision_runs.items():
        if new_date != date:
          date = new_date
          gem_plot.GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * gem_plot.GetNbinsX() + 1), '%d: %s' % (run_num, date))

    gem_plot.SetLineWidth(1)
    gem_plot.SetLineColor(1)
    gem_plot.GetXaxis().SetLabelSize(.02)
    gem_plot.GetXaxis().LabelsOption('v')

    cms_label =cms_latex()
    header = head()
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)
    gStyle.SetOptStat(0)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(gem_plot, "Muons with |#eta| > 2.5")
    #leg.Draw("same P")

    canvas.SaveAs("plots/gem_plot.pdf")

    gem_plot.Write()

    for i in range(len(run_type_plot)):

      run_type_plot[i].GetXaxis().SetTitle("Run Number")
      run_type_plot[i].GetYaxis().SetTitle("Instances")
      #run_type_plot[i].GetXaxis().SetNdivisions(5)
      
      canvas = TCanvas(run_type_plot[i].GetName() , run_type_plot[i].GetName(), 900,900)
      canvas.SetLogy()
      run_type_plot[i].Draw("HIST")
      run_type_plot[i].GetXaxis().SetLabelSize(.02)

      date = ""
      for run_num, new_date in collision_runs.items():
        if new_date != date:
          date = new_date
          run_type_plot[i].GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * run_type_plot[i].GetNbinsX() + 1), '%d: %s' % (run_num, date))
          

      run_type_plot[i].SetLineWidth(1)
      run_type_plot[i].SetLineColor(1)
      

      cms_label =cms_latex()
      header = head()
      
      gStyle.SetLegendBorderSize(0)
      gStyle.SetLegendTextSize(0.018)
      gStyle.SetOptStat(0)

      leg =TLegend(0.4,0.8,0.88,0.88);
      leg.AddEntry(run_type_plot[i], "Muons with |#eta| > 2.5")
      #leg.Draw("same P")

      canvas.SaveAs("plots/run_type_plot_%d.pdf" % i)

      run_type_plot[i].Write()

      if i == 2:
        run_type_plot[i].Divide(over_eta_int)
        run_type_plot[i].SetName('over_eta_collisions_fraction')
        run_type_plot[i].Write()

    mismatch_plot.Write()


if __name__ == "__main__":
  main()
