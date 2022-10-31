/* This script analyzed the efficiency of the EMTF NN w.r.t. displacement of the primary vertex
 *
 */

#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include "TLatex.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;

//Get distance between phi's
float DPhi(double phi1,double phi2){
  float temp=phi1-phi2;
  if (temp>3.14) temp=temp-6.28;
  if (temp<-3.14) temp=temp+6.28;
  return temp;
}


//For drawing header and title to output plots
TLatex cms_latex(){
  TLatex cms_label;
  cms_label.SetTextSize(0.04);
  cms_label.DrawLatexNDC(0.1, 0.92, "#bf{ #font[22]{CMS} #font[72]{Preliminary Simulation}}");
  return cms_label;
}

TLatex head(){
  TLatex header; 
  header.SetTextSize(0.03);
  header.DrawLatexNDC(0.63, 0.92, "#sqrt{s} = 13 TeV, Run 3 MC");
  return header; 
}


//Color set for standardized CMS plots
int DefaultColor(int j,int i){
  if (j-i==1) return 2;
  else if (j-i==2) return 4;
  else if (j-i==3) return 6;
  else if (j-i==4) return 8;
  else if (j-i==5) return 9;
  else return j;
}

//store list of root files in a given directory
void list_files(const char *dirname="C:/root/folder/", const char* ext=".root", std::vector<TString>* list = new std::vector<TString>(0)) { 
  TSystemDirectory dir(dirname, dirname); 
  TList *files = dir.GetListOfFiles(); 
  if (files) { 
    TSystemFile *file; 
    TString fname; 
    TIter next(files); 
    while ((file=(TSystemFile*)next())) { 
      fname = file->GetName(); 
      if (!file->IsDirectory() && fname.EndsWith(ext)) { 
        list->push_back((string)dirname + (string)fname.Data()); 
      } 
    } 
  } 
}

int makeEfficiencyPlotsDxy(){


  std::vector<TString> files;
  //read data
  // TString ntuple = "/eos/cms/store/user/eyigitba/emtf/matchedNtuples/matchedNtuple_HTo2LLTo4Mu_combined_cmssw_11_2_0_pre8_fwImplementation_NNv6.root";
  // TString ntuple = "/afs/cern.ch/user/e/eyigitba/L1T_dev/CMSSW_11_3_0_pre5_ptLUT/src/MatchedNtuple_muGun_v5.root";
  // TString ntuple = "/eos/cms/store/user/eyigitba/emtf/matchedNtuples/matchedNtuple_DisplacedMuGun_flatPt2to1000_posEndcap_flatXYZEtaPhi_11_3_0_pre5_10M.root";
  // TString ntuple = "/eos/cms/store/user/eyigitba/emtf/matchedNtuples/matchedNtuple_DisplacedMuGun_flatPt2to1000_allEndcap_flatXYZEtaPhi_11_3_0_pre5_NNv6_20M.root";
  // TString ntuple = "/eos/cms/store/user/eyigitba/emtf/matchedNtuples/matchedNtuple_DisplacedMuGun_flatPt2to1000_allEndcap_flatXYZEtaPhi_11_3_0_pre5_NNv7_2M.root";
  // TString ntuple = "/eos/cms/store/user/eyigitba/emtf/matchedNtuples/matchedNtuple_HTo2LongLivedTo4mu_combined_NNv6.root";
  char* dirname = "/eos/cms/store/user/eyigitba/emtf/L1Ntuples/Run3/crabOut/SingleMuon/EMTFNtuple_Run3_SingleMuon_data_13p6TeV_355872_v1/220725_163335/0000/";
  TChain * cc = new TChain("EMTFNtuple/tree");

  //add plots to event tree
  list_files(dirname, ".root", &files);

  for (int i = 0; i < files.size(); i ++)
    cc->Add(files[i]);

  TTreeReader ccReader(cc);

  //Define histograms we will be using from input files
  TTreeReaderArray<float> emtfTrackPT(ccReader,"emtfTrack_pt");
  TTreeReaderArray<float> emtfTrackPhi(ccReader,"emtfTrack_phi");
  TTreeReaderArray<float> emtfTrackEta(ccReader,"emtfTrack_eta");
  TTreeReaderArray<int> emtfTrackGMTPhi(ccReader,"emtfTrack_GMT_phi");
  TTreeReaderArray<int> emtfTrackGMTEta(ccReader,"emtfTrack_GMT_eta");
  TTreeReaderArray<short> emtfTrackSector(ccReader, "emtfTrack_sector");
  TTreeReaderArray<float> emtfTrackPTDxy(ccReader,"emtfTrack_pt_dxy");


  TTreeReaderArray<float> recoMuonPT(ccReader,"recoMuon_pt");
  TTreeReaderArray<float> recoMuonEtaSt1(ccReader,"recoMuon_etaSt1");
  TTreeReaderArray<float> recoMuonEtaSt2(ccReader,"recoMuon_etaSt2");
  TTreeReaderArray<float> recoMuonDxy(ccReader,"recoMuon_dxy");
  TTreeReaderArray<float> recoMuonPhiSt1(ccReader,"recoMuon_phiSt1");
  TTreeReaderArray<float> recoMuonPhiSt2(ccReader,"recoMuon_phiSt2");

  gStyle->SetOptStat(0);

  std::cout<<"Running on "<<cc->GetEntries()<<" evts "<<std::endl;

  //plot containers
  std::vector<TString> canvasname;
  std::vector<TString> canvasname_2D;
  std::vector<std::string> kwds;
  std::vector<TString> legs;
  std::vector<TString> legs_2D;
  std::vector<TGraphAsymmErrors*> errors;
  //std::vector<TH2F*> errors_2D;

  // cosmetic options
  std::vector<bool> grid,logY,logX;



  // initialize cuts
  float den_ptThreshold = 20.0;
  float num_ptThreshold = 10.0;
  float etaThresholdMin = 1.24;
  float etaThresholdMax = 2.5;
  float dRThreshold = .1;

  //Initialize plots
  TH1F *den_dxy = new TH1F("den_dxy", "", 20, 0, 50);
  TH1F *num_dxy = new TH1F("num_dxy", "", 20, 0, 50);
  TH1F *num_pt_dxy = new TH1F("num_dxy", "", 20, 0, 50);
  
  
  // Single muon efficiencies 
  int eventCount = 0;
  float dr; 
  float recoPhi;
  float recoEta;
  float GMTEta;
  float GMTPhi;
  float phiRad;
  int numMuons = 0;
  while(ccReader.Next()){
    eventCount++;
    if (eventCount % 1000 == 0)
      std::cout << eventCount << " events read!" << std::endl;
    

    //Find a triggered muon from event
    for(int i=0; i<recoMuonPT.GetSize(); i++){
      
      //Make sure to assign a valid coordinate to muon
      if (recoMuonPhiSt2[i] > -99)
        recoPhi = recoMuonPhiSt2[i];
      else if (recoMuonPhiSt1[i] > -99)
        recoPhi = recoMuonPhiSt1[i];
      else continue;

      if (recoMuonEtaSt2[i] > -99)
        recoEta = recoMuonEtaSt2[i];
      else if (recoMuonEtaSt1[i] > -99) 
        recoEta = recoMuonEtaSt1[i];
      else continue;
      
      //Apply trigger cuts to muon
      if (gendR[i] > dRThreshold) continue;
      if (abs(recoEta) < etaThresholdMin) continue;
      if (abs(recoEta) > etaThresholdMax) continue;
      if (recoMuonPT[i] < den_ptThreshold) continue;
      

      //Add passed to denominator w.r.t. displacement
      den_dxy->Fill(recoMuonDxy[i] * 10.0);
      numMuons ++;

      for (int j = 0; j < emtfTrackPT.GetSize(); j++) {
        
        

        //match reco-muon to a tag
        phiRad = emtfTrackPhi[j] * 3.1415 / 180.0;
        dr = TMath::Sqrt((recoEta-emtfTrackEta[j])*(recoEta-emtfTrackEta[j])+DPhi(recoPhi,phiRad)*DPhi(recoPhi,phiRad));
        if (dr > dRThreshold) continue;

        //Apply PT cut to track, see if it was triggered upon
        if (emtfTrackPT[j] > num_ptThreshold || emtfTrackPTDxy[j] > num_ptThreshold) {

          //get estimate of NN's performance w.r.t. displacement
          if (emtfTrackPTDxy[j] > num_ptThreshold)
            num_pt_dxy->Fill(recoMuonDxy[i]*10.0);

          //get estimate of BDT's performance w.r.t. displacement
          if (emtfTrackPT[j] > num_ptThreshold)
            num_dxy->Fill(recoMuonDxy[i]*10.0);
          break;
        }
        
      }
    }
  }

  // Divide histograms
  TGraphAsymmErrors* error_dxy = new TGraphAsymmErrors(num_dxy, den_dxy);
  TGraphAsymmErrors* error_pt_dxy = new TGraphAsymmErrors(num_pt_dxy, den_dxy);
  
  //Compare BDT to NN
  TString titleDxy="BDT Muon d_{xy} [cm]";
  TString titleptDxy="NN Muon d_{xy} [cm]";


  //Set labels
  error_dxy->GetXaxis()->SetTitle(titleDxy);
  error_dxy->GetYaxis()->SetTitle("L1T with BDT Efficiency");
  error_dxy->GetYaxis()->SetRangeUser(0.00001,1.2);
  error_pt_dxy->GetXaxis()->SetTitle(titleptDxy);
  error_pt_dxy->GetYaxis()->SetTitle("L1T with NN Efficiency");
  error_pt_dxy->GetYaxis()->SetRangeUser(0.00001,1.2);


  //Create legs now
  TString leg = "L1 pT > 10 GeV";
  TString leg2 = "Reco Muon > 20 GeV";
  TString leg_etaFull = "1.2 < |#eta| < 2.5";
  errors.push_back(error_dxy);
  errors.push_back(error_pt_dxy);
  
  // eta same plot
  legs.push_back(leg_etaFull);
  legs.push_back(leg);
  legs.push_back(leg2);

  legs.push_back(leg_etaFull);
  legs.push_back(leg);
  legs.push_back(leg2);

  //Create PDF output names
  canvasname.push_back("eff_dxy_pt10");
  canvasname.push_back("eff_pt_dxy_pt10");
  
  //make sure we have no repeat names
  for (int i=0; i<canvasname.size(); i++){

    //create canvas and save histos
    std::cout << "Drawing: " << errors[i]->GetName() << std::endl;
    TCanvas* c1 = new TCanvas(canvasname[i],canvasname[i],700,700);
    errors[i]->Draw("A P");
    errors[i]->SetLineWidth(3);
    errors[i]->SetLineColor(1);
    TLatex cms_label=cms_latex();
    TLatex header=head();

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(0.018); 
    TLegend * leg =new TLegend(0.4,0.8,0.88,0.88);
    leg->AddEntry(errors[i],legs[i]);

    //put histos with same cnvas name in same plot
    for (int j=i+1; j<errors.size(); j++){
      if (canvasname[i]!=canvasname[j]) continue;
      errors[j]->Draw("sames P");
      canvasname[j]="Canvas_name_already_used_action_skipping";
      // printLeg=true;
      errors[j]->SetLineWidth(3);
      errors[j]->SetLineColor(2);
      errors[j]->SetLineColor(DefaultColor(j,i));
      leg->AddEntry(errors[j],legs[j]);
    }
    leg->Draw("sames");

    c1->SaveAs(canvasname[i]+"_muGunFlatXYZ_denomModeCut.pdf");
    canvasname[i]="Canvas_name_already_used_action_skipping";

  }
  std::cout << "There are " << numMuons << " muons" << std::endl;
  return 0;
 }
