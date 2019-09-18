#define MuonTnP_cxx
#include "MuonTnP.h"
#include<iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1D.h"
#include "TLorentzVector.h"

void MuonTnP::Loop(float tagPt, float probePt, float mllDiff, float tagIso, float probeIso)
{
//   In a ROOT session, you can do:
//      root> .L MuonTnP.C
//      root> MuonTnP t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   float ptBins[] = {20, 25, 30, 40, 50, 60, 120, 200};
   float etaBins[] = {-2.4, -2.3, -2.2, -2.1, -2.0, -1.7, -1.6, -1.5, -1.4, -1.2, -0.8, -0.5, -0.3, -0.2, 0.0, 0.2, 0.3, 0.5, 0.8, 1.2, 1.4,
                      1.5, 1.6, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4};

   TH1D* hmuPt_probe_den  = new TH1D("hmuPt_probe_den", "Probe muon pt;p_{T};", 7, ptBins);
   TH1D* hmuEta_probe_den  = new TH1D("hmuEta_probe_den", "Probe muon eta;#eta;", 28, etaBins);
   TH1D* htnpMass_den  = new TH1D("htnpMass_den", "TnP pair mass", 60, 60., 120.);

   TH1D* hmuPt_probe_num_looseId  = new TH1D("hmuPt_probe_num_looseId", "Probe muon pt;p_{T};", 7, ptBins);
   TH1D* hmuEta_probe_num_looseId  = new TH1D("hmuEta_probe_num_looseId", "Probe muon eta;#eta;", 28, etaBins);
   TH1D* htnpMass_num_looseId  = new TH1D("htnpMass_num_looseId", "TnP pair mass", 60, 60., 120.);

   TH1D* hmuPt_probe_num_mediumId  = new TH1D("hmuPt_probe_num_mediumId", "Probe muon pt;p_{T};", 7, ptBins);
   TH1D* hmuEta_probe_num_mediumId  = new TH1D("hmuEta_probe_num_mediumId", "Probe muon eta;#eta;", 28, etaBins);
   TH1D* htnpMass_num_mediumId  = new TH1D("htnpMass_num_mediumId", "TnP pair mass", 60, 60., 120.);

   TH1D* hmuPt_probe_num_tightId  = new TH1D("hmuPt_probe_num_tightId", "Probe muon pt;p_{T};", 7, ptBins);
   TH1D* hmuEta_probe_num_tightId  = new TH1D("hmuEta_probe_num_tightId", "Probe muon eta;#eta;", 28, etaBins);
   TH1D* htnpMass_num_tightId  = new TH1D("htnpMass_num_tightId", "TnP pair mass", 60, 60., 120.);
  
   TH1D* hmuPt_probe_num_iso  = new TH1D("hmuPt_probe_num_iso", "Probe muon pt;p_{T};", 7, ptBins);
   TH1D* hmuEta_probe_num_iso  = new TH1D("hmuEta_probe_num_iso", "Probe muon eta;#eta;", 28, etaBins);
   TH1D* htnpMass_num_iso  = new TH1D("htnpMass_num_iso", "TnP pair mass", 60, 60., 120.);



   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   unsigned int maxEntries = nentries;
   //maxEntries = 100;

   for (Long64_t jentry=0; jentry<maxEntries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry % 10000 == 0)  std::cout << "Processed entries:" << jentry << std::endl; 
      //Loop over muons
      for (unsigned int i = 0; i < nmuon; i++) {

        //Let's look for the tag muon
        if( muon_pt->at(i) <= tagPt || std::abs(muon_eta->at(i)) > 2.4 || !muon_tight->at(i)  || muon_iso->at(i) > tagIso )   continue;

        TLorentzVector tagP4;
        tagP4.SetPtEtaPhiM(muon_pt->at(i),muon_eta->at(i),muon_phi->at(i),0.105);

        for (unsigned int j = i+1; j < nmuon; j++) {
          //Look for opposite sign muons
          if(muon_charge->at(i) + muon_charge->at(j) != 0)    continue;

          //put kinematic cuts on probe
          if( muon_pt->at(j) <= probePt || std::abs(muon_eta->at(j)) > 2.4 )   continue;

          //get Lorentz vector for the probe
          TLorentzVector probeP4;
          probeP4.SetPtEtaPhiM(muon_pt->at(j),muon_eta->at(j),muon_phi->at(j),0.105);
          //compute invariantmass
          float mll = (tagP4 + probeP4).M();

          //apply mass cut
          if( std::abs(mll - 90.) > mllDiff)   continue;

          //fill denominator histograms
          hmuPt_probe_den->Fill(probeP4.Pt());    
          hmuEta_probe_den->Fill(probeP4.Eta());    
          htnpMass_den->Fill(mll);
      
      
          if( !muon_loose->at(j) )   continue;
          hmuPt_probe_num_looseId->Fill(probeP4.Pt());
          hmuEta_probe_num_looseId->Fill(probeP4.Eta());
          htnpMass_num_looseId->Fill(mll);

          if( !muon_medium->at(j) )   continue;
          hmuPt_probe_num_mediumId->Fill(probeP4.Pt());
          hmuEta_probe_num_mediumId->Fill(probeP4.Eta());
          htnpMass_num_mediumId->Fill(mll);

          if( !muon_tight->at(j) )   continue;
          hmuPt_probe_num_tightId->Fill(probeP4.Pt());
          hmuEta_probe_num_tightId->Fill(probeP4.Eta());
          htnpMass_num_tightId->Fill(mll);

          if( muon_iso->at(j) >= probeIso )    continue;
          hmuPt_probe_num_iso->Fill(probeP4.Pt());
          hmuEta_probe_num_iso->Fill(probeP4.Eta());
          htnpMass_num_iso->Fill(mll);
  
        }//end loop over probes
      }//end loop over tags
   }//end event loop


  //Save histograms to file
  TFile* fout = TFile::Open("tnpHistos.root", "recreate"); 
  fout->cd();
  hmuPt_probe_den->Write();
  hmuEta_probe_den->Write();
  htnpMass_den->Write();
  
  hmuPt_probe_num_looseId->Write();
  hmuEta_probe_num_looseId->Write();
  htnpMass_num_looseId->Write();
  
  hmuPt_probe_num_mediumId->Write();
  hmuEta_probe_num_mediumId->Write();
  htnpMass_num_mediumId->Write();
  
  hmuPt_probe_num_tightId->Write();
  hmuEta_probe_num_tightId->Write();
  htnpMass_num_tightId->Write();
  
  hmuPt_probe_num_iso->Write();
  hmuEta_probe_num_iso->Write();
  htnpMass_num_iso->Write();
  fout->Save();
  fout->Close();
 
}
