// -*- C++ -*-
//
// Package:    MuonPOGLongExercise/MuonTnPMCtruth
// Class:      MuonTnPMCtruth
//
/**\class MuonTnPMCtruth MuonTnPMCtruth.cc MuonPOGLongExercise/MuonTnPMCtruth/plugins/MuonTnPMCtruth.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Suvankar Roy Chowdhury
//         Created:  Mon, 29 Jul 2019 14:07:38 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"
#include<vector>
#include "TGraphAsymmErrors.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class MuonTnPMCtruth : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonTnPMCtruth(const edm::ParameterSet&);
      ~MuonTnPMCtruth();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      double getmuRelIso(const pat::Muon& muon);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;  
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
     //histos to book
     TH1D* hmuPt_den;
     TH1D* hmuEta_den;
     TH1D* hmuPhi_den;
     TH1D* hnVtx_den;

     TH1D* hmuPt_num_glbMu;
     TH1D* hmuEta_num_glbMu;
     TH1D* hmuPhi_num_glbMu;
     TH1D* hnVtx_num_glbMu;
     TGraphAsymmErrors *gae_rec_pt_glbMu; 

     TH1D* hmuPt_num_looseId;
     TH1D* hmuEta_num_looseId;
     TH1D* hmuPhi_num_looseId;
     TH1D* hnVtx_num_looseId;
     TGraphAsymmErrors *gae_rec_pt_looseId; 

     TH1D* hmuPt_num_mediumId;
     TH1D* hmuEta_num_mediumId;
     TH1D* hmuPhi_num_mediumId;
     TH1D* hnVtx_num_mediumId;
     TGraphAsymmErrors *gae_rec_pt_mediumId; 

     TH1D* hmuPt_num_tightId;
     TH1D* hmuEta_num_tightId;
     TH1D* hmuPhi_num_tightId;
     TH1D* hnVtx_num_tightId;
     TGraphAsymmErrors *gae_rec_pt_tightId; 

     TH1D* hmuPt_num_iso;
     TH1D* hmuEta_num_iso;
     TH1D* hmuPhi_num_iso;
     TH1D* hnVtx_num_iso;
     TGraphAsymmErrors *gae_rec_pt_iso; 

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonTnPMCtruth::MuonTnPMCtruth(const edm::ParameterSet& iConfig)
 :
 muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
 genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
 vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc")))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  hmuPt_den  = fs->make<TH1D>("hmuPt_den", "Gen muon pt", 100, 0., 100.);
  hmuEta_den = fs->make<TH1D>("hmuEta_den", "Gen muon Eta", 25, -2.5, 2.5);
  hmuPhi_den = fs->make<TH1D>("hmuPhi_den", "Gen muon Phi", 24, -TMath::Pi(), TMath::Pi());
  hnVtx_den  = fs->make<TH1D>("hnVtx_den", "Number of good vertices", 100, 0, 100);

  hmuPt_num_glbMu  = fs->make<TH1D>("hmuPt_num_glbMu", "Matched gen muon pt", 100, 0., 100.);
  hmuEta_num_glbMu = fs->make<TH1D>("hmuEta_num_glbMu", "Matched gen muon Eta", 25, -2.5, 2.5);
  hmuPhi_num_glbMu = fs->make<TH1D>("hmuPhi_num_glbMu", "Matched gen muon Phi", 24, -TMath::Pi(), TMath::Pi());
  hnVtx_num_glbMu  = fs->make<TH1D>("hnVtx_num_glbMu", "Number of good vertices(gen-reco match)", 100, 0, 100);
  gae_rec_pt_glbMu = fs->make<TGraphAsymmErrors>();

  hmuPt_num_looseId  = fs->make<TH1D>("hmuPt_num_looseID", "Matched gen muon pt", 100, 0., 100.);
  hmuEta_num_looseId = fs->make<TH1D>("hmuEta_num_looseID", "Matched gen muon Eta", 25, -2.5, 2.5);
  hmuPhi_num_looseId = fs->make<TH1D>("hmuPhi_num_looseID", "Matched gen muon Phi", 24, -TMath::Pi(), TMath::Pi()); 
  hnVtx_num_looseId  = fs->make<TH1D>("hnVtx_num_looseID", "Number of good vertices(gen-reco match)", 100, 0, 100); 
  gae_rec_pt_looseId = fs->make<TGraphAsymmErrors>();

  hmuPt_num_mediumId  = fs->make<TH1D>("hmuPt_num_mediumId", "Matched gen muon pt", 100, 0., 100.);
  hmuEta_num_mediumId = fs->make<TH1D>("hmuEta_num_mediumId", "Matched gen muon Eta", 25, -2.5, 2.5);
  hmuPhi_num_mediumId = fs->make<TH1D>("hmuPhi_num_mediumId", "Matched gen muon Phi", 24, -TMath::Pi(), TMath::Pi());
  hnVtx_num_mediumId  = fs->make<TH1D>("hnVtx_num_mediumId", "Number of good vertices(gen-reco match)", 100, 0, 100);
  gae_rec_pt_mediumId = fs->make<TGraphAsymmErrors>();

  hmuPt_num_tightId  = fs->make<TH1D>("hmuPt_num_tightId", "Matched gen muon pt", 100, 0., 100.);
  hmuEta_num_tightId = fs->make<TH1D>("hmuEta_num_tightId", "Matched gen muon Eta", 25, -2.5, 2.5);
  hmuPhi_num_tightId = fs->make<TH1D>("hmuPhi_num_tightId", "Matched gen muon Phi", 24, -TMath::Pi(), TMath::Pi());
  hnVtx_num_tightId  = fs->make<TH1D>("hnVtx_num_tightId", "Number of good vertices(gen-reco match)", 100, 0, 100);
  gae_rec_pt_tightId = fs->make<TGraphAsymmErrors>();

  hmuPt_num_iso  = fs->make<TH1D>("hmuPt_num_iso", "Matched gen muon pt", 100, 0., 100.);
  hmuEta_num_iso = fs->make<TH1D>("hmuEta_num_iso", "Matched gen muon Eta", 25, -2.5, 2.5);
  hmuPhi_num_iso = fs->make<TH1D>("hmuPhi_num_iso", "Matched gen muon Phi", 24, -TMath::Pi(), TMath::Pi());
  hnVtx_num_iso  = fs->make<TH1D>("hnVtx_num_iso", "Number of good vertices(gen-reco match)", 100, 0, 100);
  gae_rec_pt_iso = fs->make<TGraphAsymmErrors>();

}


MuonTnPMCtruth::~MuonTnPMCtruth()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonTnPMCtruth::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByToken(vertexToken_, primaryVertices);
  //loop over vertices
  int nGoodVtx = 0;
  int selVtxindex = -1;
  for(unsigned int i = 0; i < primaryVertices->size() ;  i++) {
    const reco::Vertex& vertex = primaryVertices->at(i);
    if (!vertex.isFake() && vertex.ndof() > 4 && 
        vertex.position().Rho()<2 && std::abs(vertex.position().Z())<24. ) {
      nGoodVtx++;
      if(selVtxindex < 0) selVtxindex = i;
    } 
  }
  if(nGoodVtx == 0)   return;
  const reco::Vertex& selectedGoodVtx = primaryVertices->at(selVtxindex);

  //Loop over genparticles
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);

  std::vector<reco::GenParticle>   selectedgenMuons;
  for(auto&& genPart : *(genParticles.product())) { 
    // Check if it's a muon from Drell-Yan process
    if(genPart.isPromptFinalState() && 
       std::abs(genPart.pdgId()) == 13 && 
       genPart.status() == 1) {
      // Only muons within acceptance and pt>20 
      if(genPart.pt()>20. && std::abs(genPart.eta())<2.4) {
        hmuPt_den->Fill(genPart.pt()); 
        hmuEta_den->Fill(genPart.eta()); 
        hmuPhi_den->Fill(genPart.phi()); 
        hnVtx_den->Fill(nGoodVtx);
        selectedgenMuons.emplace_back(genPart);
      }
    }
  } 

  //Loop over reconstructed muons
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  int matchedMu= 0;
  for (auto&& muon : *(muons.product())) {
    int idx(-1), bestidx(-1);
    float bestdr(9999.);
    for (auto&& genmu : selectedgenMuons) {
      idx++; 
      float dr = deltaR(muon, genmu);
      if (dr < 0.3 && dr < bestdr) {
         bestidx = idx;
         bestdr = dr;
      }      
    }
    //std::cout << "Best idx=" << bestidx << std::endl;
    if(bestidx == -1)    continue;
    //std::cout << "Filling Histos\n";

    //Fill numerator histos for global muon
    if(!muon.isGlobalMuon())    continue;
    matchedMu++;
    hmuPt_num_glbMu->Fill(selectedgenMuons[bestidx].pt()); 
    hmuEta_num_glbMu->Fill(selectedgenMuons[bestidx].eta()); 
    hmuPhi_num_glbMu->Fill(selectedgenMuons[bestidx].phi()); 
    hnVtx_num_glbMu->Fill(nGoodVtx);

    //Fill histos for id selection
    if( !muon.isLooseMuon() )   continue;
    hmuPt_num_looseId->Fill(selectedgenMuons[bestidx].pt());
    hmuEta_num_looseId->Fill(selectedgenMuons[bestidx].eta());
    hmuPhi_num_looseId->Fill(selectedgenMuons[bestidx].phi());
    hnVtx_num_looseId->Fill(nGoodVtx); 

    if( !muon.isMediumMuon() )   continue;
    hmuPt_num_mediumId->Fill(selectedgenMuons[bestidx].pt());
    hmuEta_num_mediumId->Fill(selectedgenMuons[bestidx].eta());
    hmuPhi_num_mediumId->Fill(selectedgenMuons[bestidx].phi());
    hnVtx_num_mediumId->Fill(nGoodVtx);
   
    if( !muon.isTightMuon(selectedGoodVtx) )   continue;
    hmuPt_num_tightId->Fill(selectedgenMuons[bestidx].pt());
    hmuEta_num_tightId->Fill(selectedgenMuons[bestidx].eta());
    hmuPhi_num_tightId->Fill(selectedgenMuons[bestidx].phi());
    hnVtx_num_tightId->Fill(nGoodVtx);

    if( getmuRelIso(muon) >= 0.2 )    continue;
    hmuPt_num_iso->Fill(selectedgenMuons[bestidx].pt());
    hmuEta_num_iso->Fill(selectedgenMuons[bestidx].eta());
    hmuPhi_num_iso->Fill(selectedgenMuons[bestidx].phi());
    hnVtx_num_iso->Fill(nGoodVtx);
  }

   //std::cout << "#Matched Muons=" << matchedMu << std::endl;
  
  	
}


// ------------ method called once each job just before starting event loop  ------------
void
MuonTnPMCtruth::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonTnPMCtruth::endJob()
{
  try{
    gae_rec_pt_glbMu->Divide(hmuPt_num_glbMu, hmuPt_den, "cl=0.683 b(1,1) mode");
    gae_rec_pt_looseId->Divide(hmuPt_num_looseId, hmuPt_den, "cl=0.683 b(1,1) mode");
    gae_rec_pt_mediumId->Divide(hmuPt_num_mediumId, hmuPt_den, "cl=0.683 b(1,1) mode");
    gae_rec_pt_tightId->Divide(hmuPt_num_tightId, hmuPt_den, "cl=0.683 b(1,1) mode");
    gae_rec_pt_iso->Divide(hmuPt_num_iso, hmuPt_den, "cl=0.683 b(1,1) mode");
  } 
  catch(cms::Exception& ex) {}
 
}

double 
MuonTnPMCtruth::getmuRelIso(const pat::Muon& muon) {
  return (muon.pfIsolationR04().sumChargedHadronPt + 
         std::max(0., muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5*muon.pfIsolationR04().sumPUPt) ) / muon.pt();
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonTnPMCtruth::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonTnPMCtruth);

