// -*- C++ -*-
//
// Package:    MuonPOGLongExercise/MuonGenMatching
// Class:      MuonGenMatching
//
/**\class MuonGenMatching MuonGenMatching.cc MuonPOGLongExercise/MuonGenMatching/plugins/MuonGenMatching.cc

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
#include "TH2.h"
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

class MuonGenMatching : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonGenMatching(const edm::ParameterSet&);
      ~MuonGenMatching();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;  
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

      //histos to book
      TH1D* hmuPt_gen;
      TH1D* hmuEta_gen;
      TH1I* hnGenmu;

      TH1D* hmuPt_gen_fid;
      TH1D* hmuEta_gen_fid;
      TH1I* hnGenmu_fid;

      TH1D* hmuPt_gen_fid_prompt;
      TH1D* hmuEta_gen_fid_prompt;
      TH1I* hnGenmu_fid_prompt;

      TH1D* hmuPt_reco;
      TH1D* hmuEta_reco;
      TH1I* hnRecomu;

      TH1D* hmuPt_reco_fid;
      TH1D* hmuEta_reco_fid;
      TH1I* hnRecomu_fid;
      TH1I* hnRecomu_genMtch;
      TH1I* hnRecomu_genMtchFrompat;

      TH2D* hnGenMuvsnRecomu;

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
MuonGenMatching::MuonGenMatching(const edm::ParameterSet& iConfig)
 :
 muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
 genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
 vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc")))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  hmuPt_gen  = fs->make<TH1D>("hmuPt_gen", "Gen muon pt", 250, 0., 500.);
  hmuEta_gen = fs->make<TH1D>("hmuEta_gen", "Gen muon Eta", 25, -2.5, 2.5);

  hmuPt_gen_fid  = fs->make<TH1D>("hmuPt_gen_fid", "Gen muon(fiducial) pt", 250, 0., 500.);
  hmuEta_gen_fid = fs->make<TH1D>("hmuEta_gen_fid", "Gen muon(fiducial) Eta", 25, -2.5, 2.5);

  hmuPt_gen_fid_prompt  = fs->make<TH1D>("hmuPt_gen_fid_prompt", "Gen muon(fiducial  + prompt) pt", 250, 0., 500.);
  hmuEta_gen_fid_prompt = fs->make<TH1D>("hmuEta_gen_fid_prompt", "Gen muon(fiducial + prompt) Eta", 25, -2.5, 2.5);
  
  hnGenmu = fs->make<TH1I>("hnGenmu", ";#Gen Muons", 10, -0.5, 9.5);
  hnGenmu_fid = fs->make<TH1I>("hnGenmu_fid", ";#Gen Muons(fiducial)", 10, -0.5, 9.5);
  hnGenmu_fid_prompt = fs->make<TH1I>("hnGenmu_fid_prompt", ";#Gen Muons(fiducial and prompt)", 10, -0.5, 9.5);

  hmuPt_reco  = fs->make<TH1D>("hmuPt_reco", "Reco muon pt", 250, 0., 500.);
  hmuEta_reco = fs->make<TH1D>("hmuEta_reco", "Reco muon Eta", 25, -2.5, 2.5);

  hmuPt_reco_fid  = fs->make<TH1D>("hmuPt_reco_fid", "Reco muon(fiducial) pt", 250, 0., 500.);
  hmuEta_reco_fid = fs->make<TH1D>("hmuEta_reco_fid", "Reco muon(fiducial) Eta", 25, -2.5, 2.5);

  hnRecomu = fs->make<TH1I>("hnRecomu", ";#Reco Muons", 10, -0.5, 9.5);
  hnRecomu_fid = fs->make<TH1I>("hnRecomu_fid", ";#Reco Muons(fiducial)", 10, -0.5, 9.5);
  hnRecomu_genMtch = fs->make<TH1I>("hnRecomu_genMatch", ";#Reco Muons which are gen matched", 10, -0.5, 9.5);
  hnRecomu_genMtchFrompat = fs->make<TH1I>("hnRecomu_genMatchFrompat", ";#Reco Muons which are gen matched", 10, -0.5, 9.5);

  hnGenMuvsnRecomu = fs->make<TH2D>("hnGenMuvsnRecomu", ";nGenMuons;nRecoMuons matched to Gen", 10, -0.5, 9.5, 10, -0.5, 9.5);
}


MuonGenMatching::~MuonGenMatching()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
MuonGenMatching::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByToken(vertexToken_, primaryVertices);
  //if(nGoodVtx == 0)   return;
  const reco::Vertex& selectedGoodVtx = primaryVertices->front();

  //Loop over genparticles and plot pt/eta of gen muons with status ==1 
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);
  std::vector<reco::GenParticle> stableGenMuons; 
  unsigned int nGenmu = 0;
  unsigned int nGenmu_fid = 0;
  unsigned int nGenmu_fid_prompt = 0;

  for(auto&& genPart : *(genParticles.product())) { 
    if(std::abs(genPart.pdgId()) == 13 && genPart.status() == 1) {
      hmuPt_gen->Fill(genPart.pt()); 
      hmuEta_gen->Fill(genPart.eta());
      nGenmu++; 
      stableGenMuons.emplace_back(genPart);
      if(genPart.pt() > 20. && std::abs(genPart.eta()) < 2.4) {
        hmuPt_gen_fid->Fill(genPart.pt()); 
        hmuEta_gen_fid->Fill(genPart.eta());
        nGenmu_fid++;
        if(genPart.isPromptFinalState()) {
          nGenmu_fid_prompt++;
          hmuPt_gen_fid_prompt->Fill(genPart.pt()); 
          hmuEta_gen_fid_prompt->Fill(genPart.eta());
        }
      }
    }
  } 

  hnGenmu->Fill(nGenmu);
  hnGenmu_fid->Fill(nGenmu_fid);
  hnGenmu_fid_prompt->Fill(nGenmu_fid_prompt);
  
  //Loop over reconstructed muons to fill the properties
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  unsigned int nGenmatchedMuons = 0;
  unsigned int nGenMuonFrompat = 0;
  unsigned int nfiducialMuons = 0;
  for (auto&& muon : *(muons.product()) ) {
    hmuPt_reco->Fill(muon.pt()); 
    hmuEta_reco->Fill(muon.eta()); 
    float bestdr(9999.);    
    //loop over stable gen muons
    for(auto& genMu : stableGenMuons) {
      float dr = deltaR(muon, genMu);
      if(dr < bestdr)  bestdr = dr;
    }
    if(bestdr < 0.3)  nGenmatchedMuons++;
    //now patmuons have gen particle info embedded in them.
    const reco::GenParticle* gp = muon.genParticle();
    if(gp!=0) {
      if(std::abs(gp->pdgId()) == 13 && gp->status() == 1)   nGenMuonFrompat++;
    }
    if(muon.pt() > 20. && std::abs(muon.eta()) < 2.4) {
      hmuPt_reco_fid->Fill(muon.pt()); 
      hmuEta_reco_fid->Fill(muon.eta());
      nfiducialMuons++;
    } 
  }  

  hnRecomu_genMtch->Fill(nGenmatchedMuons);
  hnRecomu_genMtchFrompat->Fill(nGenMuonFrompat);
  hnRecomu->Fill(muons->size());
  hnRecomu_fid->Fill(nfiducialMuons);
  hnGenMuvsnRecomu->Fill(nGenmu, nGenmatchedMuons);	
}


// ------------ method called once each job just before starting event loop  ------------
void
MuonGenMatching::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonGenMatching::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonGenMatching::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonGenMatching);

