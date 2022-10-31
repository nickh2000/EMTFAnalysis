import os
from requests_futures.sessions import FuturesSession
import requests
import io
import cgi
from ROOT import *
from bs4 import BeautifulSoup
import json
import runregistry



'''This script is used to find patterns in DQM over time; Uses JSONs retrieved from DQM Online page'''

##########Define Main Builk of Code#######################
def main():

    ##Configure session with appropriate SSL certificates
    sess = FuturesSession()
    sess.verify = "/home/nick/.globus/CERN_Root_CA.crt"
    sess.cache = "/home/nick/.globus"
    sess.cert = ("/home/nick/.globus/usercert.pem", "/home/nick/.globus/userkey.pem")

    TIMEOUT = 5

    def _fetch_dqm_rows(url, timeout=TIMEOUT):
            """Return a future of DQMRows of a DQM page at url.
            Access the array of DQMRows at _resolve(self._fetch_dqm_rows(...)).data"""

            return sess.get(url, timeout=timeout, verify = "/home/nick/.globus/CERN_Root_CA.crt", stream=True)
    

    ##Define Runs of Interest######################################
    #min_run = 355100
    #min_run = 356169 # July 25
    # max_run = 357329 # August 11

    #min_run = 357330 #August 12
    #min_run = 357910 #August 23, LHC shutfown
    #max_run = 359763 #October 4

    min_run = 360129
    max_run = 360988

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
    more_high_eta = TH1D('more_high_eta', '', 100, min_run, max_run + 1)
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
    
    #Get list of runs in range of type cosmic or Commissioning
    request = runregistry.get_runs(filter={
    'class': {'or': ['Commissioning22', 'Cosmics22']},
    'run_number':{
      'and':[
        {'>=': min_run},
        {'<=': max_run},
      ]
    }
  })

    cosmics_runs = {run['run_number']: run['oms_attributes']['end_time'][5:10] for run in request}


    #Store interesting runs in a text file
    run_list = open("./plots/run_list.txt", 'w')

    #Go through all runs in range
    for run in range(min_run, max_run):
        #Print to show progress
        if not run % 100:
            print('Processing run %d' % run)

        #Get hisogram json for uGMT eta's
        eta_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1T/L1TStage2uGMT/ugmtMuonEta?formatted=true' % run
        

        response = _fetch_dqm_rows(eta_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})
        
        data = json.loads(elem.text)


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

          #if instances_neg + instances_pos > 0: print(run)
          over_eta_int.Fill(run, data['hist']['bins']['integral'])
        
        if instances_pos + instances_neg > 0: 
          print("Run Over = " + str(run))
          run_list.write(str(run) + "\n")


        #Reference json with EMTF track mode hisotra
        mode_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1T/L1TStage2EMTF/emtfTrackMode?formatted=true' % run

        response = _fetch_dqm_rows(mode_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})

        data = json.loads(elem.text)

        if data['hist'] == 'unsupported type': continue

        hist = data['hist']['bins']['content']


        #Record number of mode-zero tracks in the run, collisions runs only
        if run in collision_runs.keys(): 
          mode_zero.Fill(run, hist[0])
          if hist[0] > 0:
            run_list.write(str(run) + "\n")

        etahw_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1T/L1TStage2EMTF/MuonCand/emtfMuonhwEta?formatted=true' % run

        response = _fetch_dqm_rows(etahw_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})

        data = json.loads(elem.text)
        
        if data['hist'] == 'unsupported type': continue
        
        hist = data['hist']['bins']['content']

        #Look at number of tracks with Hardware Eta |eta| > 2.5
        if run in cosmics_runs.keys(): 
          over_etahw.Fill(run, hist[0] + hist[-1] + data['hist']['stats']['overflow'] + data['hist']['stats']['underflow'] )

        emtf_eta_path = f'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/{run}/Global/Online/ALL/L1T/L1TStage2EMTF/emtfTrackEta?formatted=true'
        response = _fetch_dqm_rows(emtf_eta_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})
        data = json.loads(elem.text)
        
        if data['hist'] == 'unsupported type': continue
        
        hist = data['hist']['bins']['content']
        

        #Store ratio of |eta| > 1.6 divided by |eta| < 1.6, Collisions only
        if run in collision_runs.keys(): 
          bound = 19
          high_eta = sum(hist[:bound]) + sum(hist[-bound:])
          low_eta = sum(hist[bound:-bound])
          if (low_eta) > 0:
            more_high_eta.Fill(run, float(high_eta) / float(low_eta))
          else:
            print("low eta was zero!")
        

        mismatch_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1TEMU/L1TdeStage2EMTF/mismatchRatio?formatted=true' % run

        response = _fetch_dqm_rows(mismatch_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})

        data = json.loads(elem.text)
        
        
        if data['hist'] == 'unsupported type': continue

        hist = data['hist']['bins']['content']


        #Find mismatch between re-emulated and hardware tracks in eta
        if run in collision_runs.keys(): 
          mismatch_plot.Fill(run, hist[2] * 100)


        
        gem_path = 'https://cmsweb.cern.ch/dqm/online/jsonfairy/archive/%d/Global/Online/ALL/L1T/L1TStage2EMTF/gemHitOccupancy?formatted=true' % run

        response = _fetch_dqm_rows(gem_path).result()
        data = str(response.text)
        soup = BeautifulSoup(data, features='lxml')
        elem = soup.find("script", {"id": "dto"})

        data = json.loads(elem.text)
        
        if data['hist'] == 'unsupported type': continue


        #Use GEM occupancy as indicator of GEM being active to compare with other errors
        if run in collision_runs.keys(): 
          gem_plot.Fill(run, data['hist']['stats']['entries'])



    #Record dates of collisions runs in plots
    date = ""
    for run_num, new_date in collision_runs.items():
      if new_date != date:
        date = new_date
        over_eta.GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * over_eta.GetNbinsX() + 1), '%d: %s' % (run_num, date))
        mode_zero.GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * mode_zero.GetNbinsX() + 1), '%d: %s' % (run_num, date))
        over_etahw.GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * over_etahw.GetNbinsX() + 1), '%d: %s' % (run_num, date))
        gem_plot.GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * gem_plot.GetNbinsX() + 1), '%d: %s' % (run_num, date))
    
    
    
    
    ############ Style Over-Eta plot#############################
    over_eta.GetXaxis().SetTitle("Run Number")
    over_eta.GetYaxis().SetTitle("Instances")
    over_eta.GetXaxis().SetLabelSize(.02)
    
    over_eta.Draw("HIST")
    over_eta.GetXaxis().SetLabelSize(.02)
    over_eta.GetXaxis().LabelsOption('v')

    over_eta.SetLineWidth(1)
    over_eta.SetLineColor(1)
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)
    gStyle.SetOptStat(0)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(over_eta, "Muons with |#eta| > 2.5")
    leg.Draw("same P")

    over_eta.Write()


    ############ Style Over-Positive-Eta plot#############################
    over_eta_pos.GetXaxis().SetTitle("Run Number")
    over_eta_pos.GetYaxis().SetTitle("Instances")
    over_eta_pos.GetXaxis().SetNdivisions(5)
    
    over_eta_pos.Draw("HIST")

    over_eta_pos.SetLineWidth(1)
    over_eta_pos.SetLineColor(1)

    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(over_eta_pos, "Muons with |#eta| > 2.5")
    leg.Draw("same P")

    over_eta_pos.Write()


    ############ Style Over-Negative-Eta plot#############################
    over_eta_neg.GetXaxis().SetTitle("Run Number")
    over_eta_neg.GetYaxis().SetTitle("Instances")
    over_eta_neg.GetXaxis().SetNdivisions(5)
    
    over_eta_neg.Draw("HIST")

    over_eta_neg.SetLineWidth(1)
    over_eta_neg.SetLineColor(1)
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(over_eta_neg, "Muons with |#eta| > 2.5")
    leg.Draw("same P")

    over_eta_neg.Write()


    ############ Style Mode-Zero plot#############################
    mode_zero.GetXaxis().SetTitle("Run Number")
    mode_zero.GetYaxis().SetTitle("Instances")
    mode_zero.GetXaxis().SetNdivisions(5)
    
    mode_zero.Draw("HIST")

    mode_zero.SetLineWidth(1)
    mode_zero.SetLineColor(1)
    mode_zero.GetXaxis().SetLabelSize(.02)
    mode_zero.GetXaxis().LabelsOption('v')

    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)
    gStyle.SetOptStat(0)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(mode_zero, "Muons with |#eta| > 2.5")
    leg.Draw("same P")

    ############ Style Over-Eta-Hardware plot#############################
    over_etahw.GetXaxis().SetTitle("Run Number")
    over_etahw.GetYaxis().SetTitle("Instances")
    over_etahw.Draw("HIST")          

    over_etahw.SetLineWidth(1)
    over_etahw.SetLineColor(1)
    over_etahw.GetXaxis().SetLabelSize(.02)
    over_etahw.GetXaxis().LabelsOption('v')
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)
    gStyle.SetOptStat(0)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(over_etahw, "Muons with |#eta| > 2.5")
    leg.Draw("same P")

    over_etahw.Write()


    ############ Style GEM-Occupancy plot#############################
    gem_plot.GetXaxis().SetTitle("Run Number")
    gem_plot.GetYaxis().SetTitle("Instances")
    gem_plot.Draw("HIST")

    gem_plot.SetLineWidth(1)
    gem_plot.SetLineColor(1)
    gem_plot.GetXaxis().SetLabelSize(.02)
    gem_plot.GetXaxis().LabelsOption('v')

    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendTextSize(0.018)
    gStyle.SetOptStat(0)

    leg =TLegend(0.4,0.8,0.88,0.88);
    leg.AddEntry(gem_plot, "Muons with |#eta| > 2.5")
    leg.Draw("same P")


    gem_plot.Write()


    #Index over-eta events by run type to compare incidence in cosmics/collisions runs
    #0 is all, 1 is cosmics, 2 is collisions
    for i in range(len(run_type_plot)):
      
      run_type_plot[i].Draw("HIST")
      run_type_plot[i].GetXaxis().SetLabelSize(.02)

      date = ""
      for run_num, new_date in collision_runs.items():
        if new_date != date:
          date = new_date
          run_type_plot[i].GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * run_type_plot[i].GetNbinsX() + 1), '%d: %s' % (run_num, date))
          

      run_type_plot[i].SetLineWidth(1)
      run_type_plot[i].SetLineColor(1)
      
      gStyle.SetLegendBorderSize(0)
      gStyle.SetLegendTextSize(0.018)
      gStyle.SetOptStat(0)

      leg =TLegend(0.4,0.8,0.88,0.88);
      leg.AddEntry(run_type_plot[i], "Muons with |#eta| > 2.5")
      leg.Draw("same P")

      run_type_plot[i].Write()


      #Normalize the collisions plot
      if i == 2:
        run_type_plot[i].Divide(over_eta_int)
        run_type_plot[i].SetName('over_eta_collisions_fraction')
        run_type_plot[i].Write()


    ############ Write Hardware-Emulator Mis-match plot#############################
    mismatch_plot.Write()




    ############ Style High-Low Eta Ratio plot#############################
    more_high_eta.GetXaxis().SetTitle("Run Number")
    more_high_eta.GetYaxis().SetTitle("Instances")
    more_high_eta.Draw("H")
    for run_num, new_date in cosmics_runs.items():
        if new_date != date:
          date = new_date
          more_high_eta.GetXaxis().SetBinLabel(int((run_num - min_run) / (max_run - min_run) * more_high_eta.GetNbinsX() + 1), '%d: %s' % (run_num, date))

    more_high_eta.Write()



#Actually run the code here    
if __name__ == "__main__":
  main()