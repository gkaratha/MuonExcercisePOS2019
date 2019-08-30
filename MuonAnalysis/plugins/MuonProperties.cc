// -*- C++ -*-
//
// Package:    MuonPOGLongExercise/MuonProperties
// Class:      MuonProperties
//
/**\class MuonProperties MuonProperties.cc MuonPOGLongExercise/MuonProperties/plugins/MuonProperties.cc

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

class MuonProperties : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonProperties(const edm::ParameterSet&);
      ~MuonProperties();
      // The elements defined inside "MuonParentage" are simply integers!
      //  - MuonParentage::PROMPT  = 0
      //  - MuonParentage::HF      = 1
      //  - MuonParentage::LF      = 2
      //  - MuonParentage::OTHER   = 3
      enum MuonParentage {PROMPT, HF, LF, OTHER};     

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      MuonParentage getParentType(const reco::GenParticle& prt);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;  
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

      bool isQCD;
     //histos to book
     TH1D* hmuPt_gen;
     TH1D* hmuEta_gen;

     TH1D* hmuPt_reco;
     TH1D* hmuEta_reco;

     TH1D* hIsGlobalMu;
     TH1D* hIsTrackerMu;
     TH1D* hglbtknormChi2;
     TH1D* hnumberOfValidMuonHits;
     TH1D* hnumberOfMatchedStations;
     TH1D* hnumberOfValidPixelHits;
     TH1D* hntrackerLayersWithMeasurement;
     TH1D* hdxy;
     TH1D* hdz;
     TH1D* hIso_sumChargedHadronPt;
     TH1D* hIso_sumNeutralHadronEt;
     TH1D* hIso_sumPhotonEt;     
     TH1D* hIso_sumPUPt;     
     TH1D* hcombinedRelIso04;
     TH1D* hdbCorrcombinedRelIso04;

     TH1I* hmuonParentage;     
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
MuonProperties::MuonProperties(const edm::ParameterSet& iConfig)
 :
 muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
 genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
 vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
 isQCD(iConfig.getUntrackedParameter<bool>("isQCD"))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  hmuPt_gen  = fs->make<TH1D>("hmuPt_gen", "Gen muon pt", 100, 0., 100.);
  hmuEta_gen = fs->make<TH1D>("hmuEta_gen", "Gen muon Eta", 25, -2.5, 2.5);

  hmuPt_reco  = fs->make<TH1D>("hmuPt_reco", "Reco muon pt", 100, 0., 100.);
  hmuEta_reco = fs->make<TH1D>("hmuEta_reco", "Reco muon Eta", 25, -2.5, 2.5);

  hIsGlobalMu = fs->make<TH1D>("hIsGlobalMu", "Is Global Muon", 2, -0.5, 1.5);
  hIsTrackerMu = fs->make<TH1D>("hIsTrackerMu", "Is Tracker Muon", 2, -0.5, 1.5);
  hglbtknormChi2 = fs->make<TH1D>("hglbtknormChi2", "Glb-track #chi^{2};Events", 55,  0.  ,  11.  );
  hnumberOfValidMuonHits = fs->make<TH1D>("hnumberOfValidMuonHits", ";Muon hits;Events", 40, -0.5 ,  39.5);
  hnumberOfMatchedStations = fs->make<TH1D>("hnumberOfMatchedStations", ";Muon stations;Events",  6, -0.5 ,   5.5);
  hnumberOfValidPixelHits = fs->make<TH1D>("hnumberOfValidPixelHits", ";Pixel hits;Events",  6, -0.5 ,   5.5);
  hntrackerLayersWithMeasurement = fs->make<TH1D>("hntrackerLayersWithMeasurement", ";Tracker layers;Events", 25, -0.5 ,  24.5 );
  hdxy = fs->make<TH1D>("hdxy", ";Track d_{xy} [cm];Events", 40, -0.10,   0.11);
  hdz = fs->make<TH1D>("hdz", ";Track d_{xy} [z];Events", 96, -24.,   24);
  hIso_sumChargedHadronPt = fs->make<TH1D>("hIso_sumChargedHadronPt", ";Charged hadron isolation;Events", 20,  0.  ,  10.  );
  hIso_sumNeutralHadronEt = fs->make<TH1D>("hIso_sumNeutralHadronEt", ";Neutral hadron isolation;Events", 20,  0.  ,  10.  );
  hIso_sumPhotonEt = fs->make<TH1D>("hIso_sumPhotonEt", ";Photon isolation ;Events", 20,  0.  ,  10.  );     
  hIso_sumPUPt = fs->make<TH1D>("hIso_sumPUPt", ";PU charged hadron isolation;Events", 20,  0.  ,  10. );     
  hcombinedRelIso04 = fs->make<TH1D>("hcombinedRelIso04", ";Combined relative isolation;Events", 100,  0.  ,   1 );     
  hdbCorrcombinedRelIso04 = fs->make<TH1D>("hdbCorrectedcombinedRelIso04", ";#Delta#beta corrected Combined relative isolation;Events", 100,  0.  ,   1. );     
  hmuonParentage = fs->make<TH1I>("hmuonParentage", "Muon parentage(PROMPT, HF, LF, OTHER);Muon Parentage;Events", 4,  -0.5, 3.5 );
}


MuonProperties::~MuonProperties()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method to determine whether a nonprompt muon is from HF  ------------
MuonProperties::MuonParentage MuonProperties::getParentType (const reco::GenParticle& prt)
{
  // Prompt: mother ID = 13, 15 (with prompt tau), 23, 24, 25 
  if(prt.isPromptFinalState() || prt.isDirectPromptTauDecayProductFinalState()) {
    return MuonParentage::PROMPT;
  }

  auto mom = &prt;
  bool sameid = true; 
  while(mom->numberOfMothers()>0 && sameid) {
    for(size_t im=0; im<mom->numberOfMothers(); ++im) {
      mom = dynamic_cast<const reco::GenParticle*>(mom->mother(im));
      if(mom->pdgId()!=prt.pdgId() && std::abs(mom->pdgId())!=15) { // if it's a tau, it's not a prompt tau -> check tau's mother
	sameid = false; 
	break;
      }
    }
  }

  // Classification
  unsigned int pdgId = abs(mom->pdgId());

  // - Prompt -- send warning: why didn't it pass isPromptFinalState() or isDirectPromptTauDecayProductFinalState()??
  if (pdgId == 13 || pdgId == 15 || pdgId == 23 || pdgId == 24 || pdgId == 25) {
    std::cout << "WARNING: mother's PDG ID is "
	      << pdgId << ", yet daughter did not pass isPromptFinalState() nor "
	      << "isDirectPromptTauDecayProductFinalState()" << std::endl;
    return MuonParentage::PROMPT;
  }

  // - From heavy-flavor hadron decay
  bool ishf = false;
  int thirdDigit = (pdgId/100) % 10;
  if ( thirdDigit == 5 || thirdDigit == 4 ) ishf = true ; // should catch all B and C mesons
  int fourthDigit = (pdgId/1000) % 10;
  if ( fourthDigit == 5 || fourthDigit == 4 ) ishf = true ; // should catch all B and C baryons
  if (ishf) return MuonParentage::HF;

  // - From light-flavor hadron decay
  return MuonParentage::LF;
}

//
// member functions
//

// ------------ method called for each event  ------------
void
MuonProperties::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  //Loop over genparticles and plot pt/eta of gen muons with status ==1 
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);

  for(auto&& genPart : *(genParticles.product())) { 
    // Check if it's a muon from Drell-Yan process
    if(std::abs(genPart.pdgId()) == 13 && genPart.status() == 1) {
      hmuPt_gen->Fill(genPart.pt()); 
      hmuEta_gen->Fill(genPart.eta()); 
    }
  } 

  //Loop over reconstructed muons to fill the properties
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  for (auto&& muon : *(muons.product()) ) {
    // Require muon to have a silicon track (i.e. global- or tracker-muon) 
    if (!muon.innerTrack().isNonnull()) continue;

    // Require that muon be within the tracker volume and have pt > 20 GeV
    if (muon.pt() <= 20. || std::abs(muon.eta()) > 2.5)   continue;
    
       // Let's check the origin of the muon: prompt, HF, LF, other?
    MuonParentage parentage = MuonParentage::OTHER;

    // For this, we need to find the gen particle associated with our reconstructed muon
    // Since we are using pat::Muons, the gen-matching is already done for us!
    // No need to loop over the GenParticle collection and perform a geometrical matching 
    const reco::GenParticle* gp = muon.genParticle();

    // Check if the pat::Muon has a GenParticle associated
    // In what cases there is no gen-matching? 
    if(gp!=0) {
      // This function determines the muon origin for you
      // Take some time to understand what it does 
      parentage = getParentType(*gp); 
    }
    else {
      // If there is no genParticle() and it's a Drell-Yan or top-top sample, stop here!
      // Proceed with the classification only when running on the QCD sample! 
      // In all the other samples, muons from light flavor decays are NOT saved
      // in the GenParticle collection. Therefore, light-flavor decays would be
      // classified as "other", which is not correct.
      // QCD samples, on the other hand, have gen-particles also for light-flavor decays 
      if(!isQCD) continue; 
    }
    hmuonParentage->Fill(parentage);

    //fill muon properties
     hmuPt_reco->Fill(muon.pt());
     hmuEta_reco->Fill(muon.eta());

     hIsGlobalMu->Fill(muon.isGlobalMuon());
     hIsTrackerMu->Fill(muon.isTrackerMuon());
     if( muon.isGlobalMuon() && muon.globalTrack().isNonnull() ) {
       hglbtknormChi2->Fill( std::min(muon.globalTrack()->normalizedChi2(), 10.99) );
       hnumberOfValidMuonHits->Fill( std::min(muon.globalTrack()->hitPattern().numberOfValidMuonHits(), 39) );
     }
     hnumberOfMatchedStations->Fill( std::min(muon.numberOfMatchedStations(), 5) );
     hnumberOfValidPixelHits->Fill( std::min(muon.innerTrack()->hitPattern().numberOfValidPixelHits(), 5) );
     hntrackerLayersWithMeasurement->Fill( std::min(muon.innerTrack()->hitPattern().trackerLayersWithMeasurement(), 14) );
     hdxy->Fill( std::min(std::max(muon.muonBestTrack()->dxy(selectedGoodVtx.position()), -0.299), 0.299) );
     hdz->Fill(std::min(std::max(muon.muonBestTrack()->dz(selectedGoodVtx.position()), -0.299), 0.299) );

     if (muon.isIsolationValid()) { 
       const reco::MuonPFIsolation &pfR04 = muon.pfIsolationR04();
       hIso_sumChargedHadronPt->Fill( std::min(pfR04.sumChargedHadronPt, (float)9.9) );
       hIso_sumNeutralHadronEt->Fill( std::min(pfR04.sumNeutralHadronEt, (float)9.9) );
       hIso_sumPhotonEt->Fill( std::min(pfR04.sumPhotonEt, (float)9.9) );
       hIso_sumPUPt->Fill( std::min(pfR04.sumPUPt, (float)9.9) );
       hcombinedRelIso04->Fill((pfR04.sumChargedHadronPt + std::max(0., double(pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt)))/muon.pt());
       hdbCorrcombinedRelIso04->Fill((pfR04.sumChargedHadronPt + std::max(0., double(pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt)))/muon.pt());
    }    

  }  
  	
}


// ------------ method called once each job just before starting event loop  ------------
void
MuonProperties::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonProperties::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonProperties::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonProperties);

