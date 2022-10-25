
// Contains all the branches in the NTuples
// Make sure it's up to date with the NTuple
#include "EMTFAnalyzer/NTupleMaker/interface/Read_FlatNtuple.h"
// #include "EMTFAnalyzer/NTupleMaker/interface/Read_FlatNtuple_Slim.h"

const bool verbose = false;
const int MAX_EVT = -1;  // Max number of events to process
const int PRT_EVT = 10000;   // Print every N events

// RECO muon eta ranges
const std::vector<float> eta_min = {1.25, 1.25, 1.55, 1.85, 2.1};
const std::vector<float> eta_max = {2.4,  1.55, 1.85, 2.1,  2.4};
// EMTF quality thresholds
const std::vector<int> qual_cuts = {12};
// EMTF pT thresholds
const std::vector<int> pt_cuts = {0, 22, 30};

// const TString in_dir = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/ntuples/HADD";
const TString in_dir = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/ntuples/SingleMuon";

// std::vector<TString> in_file_names = { "NTuple_SingleMuon_FlatNtuple_Run_2018D_v1_2018_09_18_SingleMuon_2018_emul_ph_lut_v2_coord.root",
// 					     "NTuple_SingleMuon_FlatNtuple_Run_2018D_v1_2018_09_19_SingleMuon_2018_emul_102X_ReReco_v1_321988_bugFix.root" };
std::vector<TString> in_file_names = { "FlatNtuple_Run_2018D_v1_2018_09_18_SingleMuon_2018_emul_ph_lut_v2_coord/NTuple_0.root",
				       "FlatNtuple_Run_2018D_v1_2018_09_19_SingleMuon_2018_emul_102X_ReReco_v1_321988_bugFix/NTuple_0.root" };
const TString out_dir = "plots";

const int pt_bins[3] = {50, 0, 50};

const float BIT = 0.001;  // Small shift for filling inside bin edges

TString STR(int x) { return std::to_string(x); }
TString STR(float x) { 
  std::ostringstream out;
  out << std::setprecision(2) << x;
  return out.str();
}


// Main function to make efficiency turn-on plots
void EMTF_efficiency(std::vector<TString> _in_file_names = {}, TString _label = "") {

  if (_in_file_names.size() > 0) in_file_names = _in_file_names;
   
  std::vector<TChain*> in_chains;
  // List of input files
  for (const auto & file_name : in_file_names) {
    in_chains.push_back( new TChain("FlatNtupleData/tree") );
    in_chains.back()->Add(in_dir+"/"+file_name);
  }

  std::vector<TH1F*> h_pt_num;  // Efficiency turn-on curve numerator
  std::vector<TH1F*> h_pt_den;  // Efficiency turn-on curve denominator
  std::vector<TEfficiency*> h_pt_eff;  // Efficiency turn-on curve

  // Book histograms
  for (int i = 0; i < in_chains.size(); i++) {
    for (int j = 0; j < eta_min.size(); j++) {
      for (const auto & qual_cut : qual_cuts) {
	for (const auto & pt_cut : pt_cuts) {

	  //  TString h_name  = "file_"+i+"_eta_"+int(eta_min.at(j)*100)+"_to_"+int(eta_max.at(j)*100)+"_qual_"+qual_cut+"_pt_"+pt_cut;
	  TString h_name  = "file_"+STR(i)+"_eta_"+STR(int(eta_min.at(j)*100))+"_to_"+STR(int(eta_max.at(j)*100))+"_qual_"+STR(qual_cut)+"_pt_"+STR(pt_cut);
	  TString h_title = "from file "+STR(i)+", "+STR(eta_min.at(j))+" < |#eta| < "+STR(eta_max.at(j))+", quality #geq "+STR(qual_cut)+", p_{T} #geq "+STR(pt_cut);
	  
	  h_pt_num.push_back( new TH1F( "h_pt_num_"+h_name, "RECO pT "+h_title+" (numerator)",   pt_bins[0], pt_bins[1], pt_bins[2] ) );
	  h_pt_den.push_back( new TH1F( "h_pt_den_"+h_name, "RECO pT "+h_title+" (denominator)", pt_bins[0], pt_bins[1], pt_bins[2] ) );
	  h_pt_eff.push_back( new TEfficiency( "h_pt_eff_"+h_name, "RECO pT "+h_title+" (efficiency)",  pt_bins[0], pt_bins[1], pt_bins[2] ) );

	} // End loop: for (const auto & pt_cut : pt_cuts)
      } // End loop: for (const auto & qual_cut : qual_cuts)
    } // End loop: for(int j = 0; j < eta_min.size(); j++)
  } // End loop: for (int i = 0; i < in_chains.size(); i++)

  // Number of histograms under each category
  int hChain = pt_cuts.size() * qual_cuts.size() * eta_min.size();
  int hEta   = pt_cuts.size() * qual_cuts.size();
  int hQual  = pt_cuts.size();


  // Loop over events in each chain
  for (int iChain = 0; iChain < in_chains.size(); iChain++) {

    TChain * in_chain = in_chains.at(iChain);
    int nEvents = (MAX_EVT > 0 ? MAX_EVT : in_chain->GetEntries());

    std::cout << "\n******* About to loop over " << nEvents << " events in file " << iChain << " *******" << std::endl;
    
    InitializeMaps();
    SetBranchAddresses(in_chain);
    std::cout << "Successfully set branch addresses" << std::endl;

    for (int iEvt = 0; iEvt < nEvents; iEvt++) {
      if ( (iEvt % PRT_EVT) == 0 ) {
	std::cout << "\n*************************************" << std::endl;
	std::cout << "Looking at event " << iEvt << " out of " << nEvents << std::endl;
	std::cout << "*************************************" << std::endl;
      }

      in_chain->GetEntry(iEvt);
    
      // From Read_FlatNtuple.h, use 'I("branch_name")' to get an integer branch value, 'F("branch_name") to get a float
      if (verbose) std::cout << "\n" << I("nRecoMuons") << " RECO muons in the event" << std::endl;
	
      // Skip events without HLT fired
      if ( I("nRecoMuonsTrig") == 0 ) continue;
	
      // Loop over RECO muons
      for (int iMu = 0; iMu < I("nRecoMuons"); iMu++) {   

	// Require muon to pass tight ID
	if ( I("reco_ID_tight", iMu) != 1 ) continue;
	// Require muon to be a "probe", i.e. at least one other muon fires the HLT
	if ( I("nRecoMuonsTrig") < 2 && I("reco_trig_ID", iMu) > 0 ) continue;

	for (int iEta = 0; iEta < eta_min.size(); iEta++) {
	  for (int iQual = 0; iQual < qual_cuts.size(); iQual++) {
	    for (int iPt = 0; iPt < pt_cuts.size(); iPt++) {

	      int iHist = (iChain * hChain) + (iEta * hEta) + (iQual * hQual) + iPt;

	      // Fill all RECO muons in eta range
	      if ( abs( F("reco_eta", iMu) ) > eta_min.at(iEta) && 
		   abs( F("reco_eta", iMu) ) < eta_max.at(iEta) ) {

		h_pt_den.at(iHist)->Fill( std::min( F("reco_pt", iMu), pt_bins[2] - BIT ) );

		int iTrk = I("reco_dR_match_emu_iTrk", iMu);
		if (iTrk < 0) continue;

		// Fill only RECO muons matched to EMTF track passing time, quality, and pT cuts
		if ( I("trk_BX",   iTrk) == 0 && 
		     I("trk_qual", iTrk) >= qual_cuts.at(iQual) &&
		     F("trk_pt",   iTrk) >= pt_cuts.at(iPt) ) {

		  h_pt_num.at(iHist)->Fill( std::min( F("reco_pt", iMu), pt_bins[2] - BIT ) );

		} // End if matched EMTF track passes cuts

	      } // End if RECO muon falls in eta range

	    } // End loop: for (int iPt = 0; iPt < pt_cuts.size(); iPt++)
	  } // End loop: for (int iQual = 0; iQual < qual_cuts.size(); iQual++)
	} // End loop: for (int iEta = 0; iEta < eta_min.size(); iEta++)

      } // End loop: for (int iMu = 0; iMu < I("nRecoMuons"); iMu++) {   
    } // End loop: for (int iEvt = 0; iEvt < nEvents; iEvt++) {
  } // End loop: for (const int iChain = 0; iChain < in_chains.size(); iChain++)

  std::cout << "\n******* Finished looping over events in all " << in_chains.size() << " files *******" << std::endl;

  TFile out_file(out_dir+"/EMTF_efficiency_tightID"+_label+".root", "RECREATE");
  out_file.cd();

  for (int i = 0; i < h_pt_num.size(); i++) {

    TString h_name  = h_pt_eff.at(i)->GetName();
    TString h_title = h_pt_eff.at(i)->GetTitle();
    h_pt_eff.at(i) = new TEfficiency( (*h_pt_num.at(i)), (*h_pt_den.at(i)) );
    h_pt_eff.at(i)->SetName (h_name);
    h_pt_eff.at(i)->SetTitle(h_title);
    h_pt_eff.at(i)->Write();

    h_pt_num.at(i)->SetLineWidth(2);
    h_pt_den.at(i)->SetLineWidth(2);
    h_pt_eff.at(i)->SetLineWidth(2);

    int nHists = eta_min.size() * qual_cuts.size() * pt_cuts.size();

    h_pt_num.at(i)->SetLineColor(4); // kBlue
    h_pt_den.at(i)->SetLineColor(1); // kBlack
    h_pt_eff.at(i)->SetLineColor( 2 + (i/nHists) ); // kRed, kGreen, kBlue ...

    h_pt_num.at(i)->Write();
    h_pt_den.at(i)->Write();
    h_pt_eff.at(i)->Write();

  }

  out_file.Close();
  
} // End function: void 
