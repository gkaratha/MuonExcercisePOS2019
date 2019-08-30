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
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/Common/interface/TriggerResults.h"
//#include "​DataFormats/L1Trigger/interface/Muon.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

#include "TLorentzVector.h"
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

class MuonTnP : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit MuonTnP(const edm::ParameterSet&);
  ~MuonTnP();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  double getmuRelIso(const pat::Muon& muon);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;  
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  
  edm::EDGetTokenT<edm::TriggerResults> hltToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> objectToken_;
  
  edm::EDGetTokenT<l1t::MuonBxCollection> l1MuonsToken_; 
  
  //histos to book
  TH1D* hmuPt_probe_den;
  TH1D* htnpMass_den;
  TH1D* hnVtx_den;
  
  TH1D* hmuPt_probe_num_looseId;
  TH1D* htnpMass_num_looseId;
  TH1D* hnVtx_num_looseId;
  
  TH1D* hmuPt_probe_num_mediumId;
  TH1D* htnpMass_num_mediumId;
  TH1D* hnVtx_num_mediumId;
  
  TH1D* hmuPt_probe_num_tightId;
  TH1D* htnpMass_num_tightId;
  TH1D* hnVtx_num_tightId;
  
  TH1D* hmuPt_probe_num_iso;
  TH1D* htnpMass_num_iso;
  TH1D* hnVtx_num_iso;
  
  
  std::string HLTPath_;
  
  int event,run_number,ls,nmuon,njpsi;
  float pvertex_x,pvertex_y,pvertex_z,bspot_x,bspot_y,bspot_z;
  //gen
  std::vector<float> genmuon_pt,genmuon_eta,genmuon_phi,genmuon_charge,genmuon_momId,genmuon_grandmomId;
  //trg
  bool hltpath;
  std::vector<float> l1muon_pt,l1muon_eta,l1muon_phi,hltmuon_pt,hltmuon_eta,hltmuon_phi;

  //muon
  std::vector<bool> muon_soft,muon_medium,muon_loose,muon_tight;
  std::vector<float> muon_pt,muon_eta,muon_phi,muon_charge,muon_vx,muon_vy,muon_vz,muon_iso;

  //jpsi
  std::vector<float> jpsi_pt,jpsi_eta,jpsi_phi,jpsi_ctxy,jpsi_x,jpsi_y,jpsi_z,jpsi_prob,jpsi_mass,jpsi_mass_unfit,jpsi_mu1pt,jpsi_mu1eta,jpsi_mu1phi,jpsi_mu2pt,jpsi_mu2eta,jpsi_mu2phi;

  edm::Service<TFileService> fs;
  //TTree * t1;

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
MuonTnP::MuonTnP(const edm::ParameterSet& iConfig)
 :
 muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
 vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
 hltToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltTag"))),
 objectToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjectTag"))),
 l1MuonsToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1muons")))
{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;
   //Histogramming part
   hmuPt_probe_den  = fs->make<TH1D>("hmuPt_probe_den", "Probe muon pt", 100, 0., 100.);
   htnpMass_den  = fs->make<TH1D>("htnpMass_den", "TnP pair mass", 60, 60., 120.);
   hnVtx_den  = fs->make<TH1D>("hnVtx_den", "Number of good vertices", 100, 0, 100);

   hmuPt_probe_num_looseId  = fs->make<TH1D>("hmuPt_probe_num_looseId", "Probe muon pt", 100, 0., 100.);
   htnpMass_num_looseId  = fs->make<TH1D>("htnpMass_num_looseId", "TnP pair mass", 60, 60., 120.);
   hnVtx_num_looseId  = fs->make<TH1D>("hnVtx_num_looseID", "Number of good vertices", 100, 0, 100);

   hmuPt_probe_num_mediumId  = fs->make<TH1D>("hmuPt_probe_num_mediumId", "Probe muon pt", 100, 0., 100.);
   htnpMass_num_mediumId  = fs->make<TH1D>("htnpMass_num_mediumId", "TnP pair mass", 60, 60., 120.);
   hnVtx_num_mediumId  = fs->make<TH1D>("hnVtx_num_mediumId", "Number of good vertices", 100, 0, 100);

   hmuPt_probe_num_tightId  = fs->make<TH1D>("hmuPt_probe_num_tightId", "Probe muon pt", 100, 0., 100.);
   htnpMass_num_tightId  = fs->make<TH1D>("htnpMass_num_tightId", "TnP pair mass", 60, 60., 120.);
   hnVtx_num_tightId  = fs->make<TH1D>("hnVtx_num_tightId", "Number of good vertices", 100, 0, 100);
  
   hmuPt_probe_num_iso  = fs->make<TH1D>("hmuPt_probe_num_iso", "Probe muon pt", 100, 0., 100.);
   htnpMass_num_iso  = fs->make<TH1D>("htnpMass_num_iso", "TnP pair mass", 60, 60., 120.);
   hnVtx_num_iso  = fs->make<TH1D>("hnVtx_num_iso", "Number of good vertices", 100, 0, 100);
}


MuonTnP::~MuonTnP()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonTnP::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByToken(vertexToken_, primaryVertices);
  //we can do this because we have a vertex filter in our config
  const reco::Vertex& selectedGoodVtx = primaryVertices->front();
  const int nGoodVtx = primaryVertices->size();

  //Loop over reconstructed muons
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  for (unsigned int i = 0; i < muons.product()->size(); i++) {
    const pat::Muon& tag = muons->at(i);
    if( tag.pt() <= 24. || std::abs(tag.eta()) > 2.4 || 
        !tag.isTightMuon(selectedGoodVtx)  || getmuRelIso(tag) > 0.15 )   continue;
    //once a tag is found, look for probes
    for (unsigned int j = i+1; j < muons.product()->size(); j++) {
      const pat::Muon& probe = muons->at(j);
      if(tag.charge() + probe.charge()!= 0)    continue;
      if( probe.pt() <= 20. || std::abs(tag.eta()) > 2.4 || 
          !probe.isGlobalMuon())   continue;
      float mll = (tag.p4() + probe.p4()).M();
      
      if( std::abs(mll - 90.) > 10.)   continue;

      hmuPt_probe_den->Fill(probe.pt());    
      htnpMass_den->Fill(mll);      
      hnVtx_den->Fill(nGoodVtx);
      
      if( !probe.isLooseMuon() )   continue;
      hmuPt_probe_num_looseId->Fill(probe.pt());
      htnpMass_num_looseId->Fill(mll);
      hnVtx_num_looseId->Fill(nGoodVtx);

      if( !probe.isMediumMuon() )   continue;
      hmuPt_probe_num_mediumId->Fill(probe.pt());
      htnpMass_num_mediumId->Fill(mll);
      hnVtx_num_mediumId->Fill(nGoodVtx);

      if( !probe.isTightMuon(selectedGoodVtx) )   continue;
      hmuPt_probe_num_tightId->Fill(probe.pt());
      htnpMass_num_tightId->Fill(mll);
      hnVtx_num_tightId->Fill(nGoodVtx);

      if( getmuRelIso(probe) >= 0.2 )    continue;
      hmuPt_probe_num_iso->Fill(probe.pt());
      htnpMass_num_iso->Fill(mll);
      hnVtx_num_iso->Fill(nGoodVtx);
      
    }
  }

  	
}


// ------------ method called once each job just before starting event loop  ------------
void
MuonTnP::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonTnP::endJob()
{
 
}

double 
MuonTnP::getmuRelIso(const pat::Muon& muon) {
  return (muon.pfIsolationR04().sumChargedHadronPt + 
         std::max(0., muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5*muon.pfIsolationR04().sumPUPt) ) / muon.pt();
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonTnP::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonTnP);

