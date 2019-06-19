// Package:    POS2019/MuonAnalysis
// Class:      MuonAnalysis
// 
/**

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <vector>
#include "TTree.h"
#include <string>
#include <iostream>
#include "TMath.h"
#include "DataFormats/Common/interface/Ref.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TLorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"



using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
template<typename T1>
class MuonAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

  typedef std::vector<T1> T1Collection;
  typedef edm::Ref<T1Collection> T1Ref;

public:
  explicit MuonAnalysis(const edm::ParameterSet&);
  ~MuonAnalysis();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void genMuon(edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > & packedGenToken_,const edm::Event& iEvent);
   

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken l1MuonsToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_; 
  edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> trigobjectsToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
  string HLTPath_;
 
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
  TTree * t1;
  bool Data=true; bool SkipEvtWithNoFire=false; double MuonPtCut=0;
  double McutMin=0; double McutMax=100; double PtMuMuCut=0;
  bool ReconstructMuMuVtx=false; bool SkipEvtWithNoVtx=false;

  
    // ----------member data ---------------------------
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
template<typename T1>
MuonAnalysis<T1>::MuonAnalysis(const edm::ParameterSet& iConfig): 
 beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot"))),
 vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
 muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
 l1MuonsToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1muons"))),
 trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag> ("triggerresults"))),
 trigobjectsToken_(consumes<vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag> ("triggerobjects"))),
 packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
HLTPath_(iConfig.getParameter<string> ("HLTPath"))
{
  edm::ParameterSet runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
  Data=runParameters.getParameter<bool>("Data");
  SkipEvtWithNoFire=runParameters.getParameter<bool>("SkipEvtWithNoFire");
  MuonPtCut=runParameters.getParameter<double>("MuonPtCut");
  McutMin=runParameters.getParameter<double>("McutMin");
  McutMax=runParameters.getParameter<double>("McutMax");
  PtMuMuCut=runParameters.getParameter<double>("PtMuMuCut");
  ReconstructMuMuVtx=runParameters.getParameter<bool>("ReconstructMuMuVtx");
  SkipEvtWithNoVtx=runParameters.getParameter<bool>("SkipEvtWithNoVtx");
}

template<typename T1>
MuonAnalysis<T1>::~MuonAnalysis()
{
  // cout<<"total "<<trg_counter<<" fires "<<fire_counter<<" l3 "<<l3_counter<<endl;
   // do anything here that needs to be done at desctruction time

}


template<typename T1>
void MuonAnalysis<T1>::genMuon(edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > & packedGenToken_,const edm::Event& iEvent){
  edm::Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_,packed);
  for ( const pat::PackedGenParticle  & gen: *packed){
    if (gen.pdgId()!=13 && gen.pdgId()!=-13) continue;
    genmuon_pt.push_back(gen.pt()); genmuon_eta.push_back(gen.eta());
    genmuon_phi.push_back(gen.phi()); genmuon_charge.push_back(gen.charge());
    genmuon_momId.push_back(gen.pdgId());
    const reco::Candidate * motherInPrunedCollection =gen.mother(0);
    if (motherInPrunedCollection!=NULL){
      genmuon_momId.push_back(motherInPrunedCollection->pdgId());
      genmuon_grandmomId.push_back(motherInPrunedCollection->mother()->pdgId());
    }
  }
}


  
template<typename T1>
void
MuonAnalysis<T1>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
  //clear vectors
  genmuon_pt.clear(),genmuon_eta.clear(),genmuon_phi.clear(),genmuon_charge.clear(),genmuon_momId.clear(),genmuon_grandmomId.clear();
  l1muon_pt.clear(),l1muon_eta.clear(),l1muon_phi.clear(),hltmuon_pt.clear(),hltmuon_eta.clear(),hltmuon_phi.clear();
  muon_soft.clear(),muon_medium.clear(),muon_loose.clear(),muon_tight.clear();
  muon_pt.clear(),muon_eta.clear(),muon_phi.clear(),muon_charge.clear(),muon_vx.clear(),muon_vy.clear(),muon_vz.clear(),muon_iso.clear();
  jpsi_pt.clear(),jpsi_eta.clear(),jpsi_phi.clear(),jpsi_ctxy.clear(),jpsi_x.clear(),jpsi_y.clear(),jpsi_z.clear(),jpsi_prob.clear(),jpsi_mass.clear(),jpsi_mass_unfit.clear(),jpsi_mu1pt.clear(),jpsi_mu1eta.clear(),jpsi_mu1phi.clear(),jpsi_mu2pt.clear(),jpsi_mu2eta.clear(),jpsi_mu2phi.clear();
  hltpath=false; njpsi=0; nmuon=0;
  //Get data
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot); 
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  //continue if there are no vertices
  if (vertices->size()==0) return;
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_,muons);
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  edm::Handle<l1t::MuonBxCollection> l1Muons;
  iEvent.getByToken(l1MuonsToken_,l1Muons);
  edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(trigobjectsToken_ ,triggerObjects);
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  //common stuff
  run_number=iEvent.id().run(); ls=iEvent.luminosityBlock();  
  bspot_x= theBeamSpot->x0(); bspot_y= theBeamSpot->y0();
  bspot_z= theBeamSpot->z0();
  reco::TrackBase::Point  vertex_point; 
  const  reco::Vertex firstGoodVertex=vertices->front();
  for (const reco::Vertex &vtx : *vertices) {
   bool isFake = vtx.isFake();
   if ( isFake || !vtx.isValid () ) continue;
   pvertex_x=vtx.x(); pvertex_y=vtx.y(); pvertex_z=vtx.z();
   vertex_point.SetCoordinates(vtx.x(),vtx.y(),vtx.z());
   break;
 }

  if (!Data) genMuon(packedGenToken_,iEvent);
  

 for(typename std::vector< l1t::Muon >::const_iterator mu=l1Muons->begin(0); mu !=l1Muons->end(0); mu++){
    l1muon_pt.push_back(mu->et()); l1muon_eta.push_back(mu->eta());
    l1muon_phi.push_back(mu->phi());
  }
 
   const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
   for(unsigned int i=0; i<triggerObjects->size(); ++i){
     pat::TriggerObjectStandAlone itrg=triggerObjects->at(i);
     itrg.unpackPathNames(trigName);
     std::vector<std::string> const& pathnames = itrg.pathNames();
     bool save=false;
     for(unsigned int ipath=0; ipath<pathnames.size(); ++ipath){
        if (pathnames[ipath].find(HLTPath_)) save=true;
     }
    if(save){
      hltpath=true;
      hltmuon_pt.push_back(itrg.pt());
      hltmuon_eta.push_back(itrg.eta());
      hltmuon_phi.push_back(itrg.phi());
    }
   }     
   
  if (SkipEvtWithNoFire && Data && !hltpath) return;
 std::vector<pat::Muon> jpsicands;
 for (const pat::Muon &mu : *muons){    
    muon_pt.push_back(mu.pt()); muon_phi.push_back(mu.phi());
    muon_eta.push_back(mu.eta()); muon_charge.push_back(mu.charge());    
    muon_vx.push_back(mu.vx()); muon_vy.push_back(mu.vy());
    muon_vz.push_back(mu.vz()); muon_medium.push_back(mu.isMediumMuon());
    muon_loose.push_back(mu.isLooseMuon());
    muon_tight.push_back(mu.isTightMuon(firstGoodVertex));
    muon_soft.push_back(mu.isSoftMuon(firstGoodVertex));
    const MuonPFIsolation&  isol=mu.pfIsolationR04();
    muon_iso.push_back((isol.sumChargedHadronPt+max(0.,isol.sumNeutralHadronEt+isol.sumPhotonEt-0.5*isol.sumPUPt))/mu.pt());
    if(mu.pt()>PtMuMuCut && mu.isSoftMuon(firstGoodVertex)){
      jpsicands.push_back(mu);  
    }
    nmuon++;
  }
   
if(jpsicands.size()>1 ) {
  for (std::vector<pat::Muon>::const_iterator mu1=jpsicands.begin(); mu1!=jpsicands.end(); ++mu1){
   for (std::vector<pat::Muon>::const_iterator mu2=mu1+1; mu2!=jpsicands.end(); ++mu2){
      if (mu1->charge()==mu2->charge()) continue;
      TLorentzVector vmu1,vmu2; 
      vmu1.SetPtEtaPhiM(mu1->pt(),mu1->eta(),mu1->phi(),0.105);
      vmu2.SetPtEtaPhiM(mu2->pt(),mu2->eta(),mu2->phi(),0.105);
      if ( (vmu1+vmu2).M()<McutMin || (vmu1+vmu2).M()>McutMax) continue;
      jpsi_mass_unfit.push_back((vmu1+vmu2).M());
      std::vector<reco::TransientTrack> jpsitrk; jpsitrk.reserve(2);
      jpsitrk.emplace_back(reco::TransientTrack(*(mu1->bestTrack()),&(*bFieldHandle)));
      jpsitrk.emplace_back(reco::TransientTrack(*(mu2->bestTrack()),&(*bFieldHandle)));
      KalmanVertexFitter vtxFitter(true);
      TransientVertex jpsivtx=vtxFitter.vertex(jpsitrk);
      if(!jpsivtx.isValid()) continue; 
      jpsi_x.push_back(jpsivtx.position().x());
      jpsi_y.push_back(jpsivtx.position().y());
      jpsi_z.push_back(jpsivtx.position().z());
      jpsi_prob.push_back(ChiSquaredProbability(jpsivtx.totalChiSquared(), jpsivtx.degreesOfFreedom()));
      std::vector<reco::TransientTrack> refited=jpsivtx.refittedTracks();
      vmu1.SetPtEtaPhiM(refited[0].track().pt(),refited[0].track().eta(),refited[0].track().phi(),0.105);
      vmu2.SetPtEtaPhiM(refited[1].track().pt(),refited[1].track().eta(),refited[1].track().phi(),0.105);
      jpsi_mu1pt.push_back(vmu1.Pt()); jpsi_mu1eta.push_back(vmu1.Eta());
      jpsi_mu1phi.push_back(vmu1.Phi()); jpsi_mu2pt.push_back(vmu2.Pt()); 
      jpsi_mu2eta.push_back(vmu2.Eta()); jpsi_mu2phi.push_back(vmu2.Phi());
      jpsi_pt.push_back((vmu1+vmu2).Pt()); jpsi_eta.push_back((vmu1+vmu2).Eta());
      jpsi_phi.push_back((vmu1+vmu2).Phi()); jpsi_mass.push_back((vmu1+vmu2).M());
      TVector2 pvxy,jpsivxy,ptj; pvxy.Set(vertex_point.x(),vertex_point.y()); 
      jpsivxy.Set(jpsivtx.position().x(),jpsivtx.position().y());  
      ptj.Set((vmu1+vmu2).Px(),(vmu1+vmu2).Py());
      float Lj=(jpsivxy-pvxy)*ptj/ptj.Mod();
      float lj=Lj*(vmu1+vmu2).M()/ptj.Mod()*10;
      jpsi_ctxy.push_back(lj);
      njpsi++;
   }
  }
 }
 if (njpsi==0 && SkipEvtWithNoVtx && Data) return;
t1->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
template<typename T1>
void 
MuonAnalysis<T1>::beginJob()
{
 t1=fs->make<TTree>("mytree","mytree");
 t1->Branch("event",&event); t1->Branch("run_number",&run_number);
 t1->Branch("ls",&ls);
 t1->Branch("pvertex_x",&pvertex_x); t1->Branch("pvertex_y",&pvertex_y);
 t1->Branch("pvertex_z",&pvertex_z); t1->Branch("bspot_x",&bspot_x);
 t1->Branch("bspot_y",&bspot_y); t1->Branch("bspot_z",&bspot_z);
 //gen
 t1->Branch("genmuon_pt",&genmuon_pt); t1->Branch("genmuon_eta",&genmuon_eta);
 t1->Branch("genmuon_phi",&genmuon_phi); t1->Branch("genmuon_charge",&genmuon_charge);
 t1->Branch("genmuon_motmId",&genmuon_momId); t1->Branch("genmuon_grandmomId",&genmuon_grandmomId);
 //trigger
 t1->Branch("hltpath",&hltpath);  
 t1->Branch("l1muon_pt",&l1muon_pt); t1->Branch("l1muon_eta",&l1muon_eta);
 t1->Branch("l1muon_phi",&l1muon_phi);
 t1->Branch("hltmuon_pt",&hltmuon_pt); t1->Branch("HLTmuon_eta",&hltmuon_eta);
 t1->Branch("hltmuon_phi",&hltmuon_phi);
 //reco 
t1->Branch("nmuon",&nmuon);
 t1->Branch("muon_pt",&muon_pt); t1->Branch("muon_eta",&muon_eta);
 t1->Branch("muon_phi",&muon_phi); t1->Branch("muon_charge",&muon_charge); 
 t1->Branch("muon_vx",&muon_vx); t1->Branch("muon_vy",&muon_vy);
 t1->Branch("muon_vz",&muon_vz); t1->Branch("muon_iso",&muon_iso);  
 t1->Branch("muon_soft",&muon_soft); t1->Branch("muon_loose",&muon_loose);
 t1->Branch("muon_medium",&muon_medium); t1->Branch("muon_tight",&muon_tight);
 //jpsi
 t1->Branch("njpsi",&njpsi);
 t1->Branch("jpsi_pt",&jpsi_pt); t1->Branch("jpsi_eta",&jpsi_eta);
 t1->Branch("jpsi_phi",&jpsi_phi);  t1->Branch("jpsi_ctxy",&jpsi_ctxy);
 t1->Branch("jpsi_x",&jpsi_x); t1->Branch("jpsi_y",&jpsi_y);
 t1->Branch("jpsi_z",&jpsi_z); t1->Branch("jpsi_prob",&jpsi_prob);
 t1->Branch("jpsi_mass",&jpsi_mass); t1->Branch("jpsi_mass_unfit",&jpsi_mass_unfit);
 t1->Branch("jpsi_mu1pt",&jpsi_mu1pt); t1->Branch("jpsi_mu1eta",&jpsi_mu1eta);
 t1->Branch("jpsi_mu1phi",&jpsi_mu1phi); t1->Branch("jpsi_mu2pt",&jpsi_mu2pt);
 t1->Branch("jpsi_mu2eta",&jpsi_mu2eta); t1->Branch("jpsi_mu2phi",&jpsi_mu2phi);
}

// ------------ method called once each job just after ending the event loop  ------------
template<typename T1>
void 
MuonAnalysis<T1>::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<typename T1>
void
MuonAnalysis<T1>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


///////////////////////
  
//define this as a plug-in
typedef MuonAnalysis<reco::RecoCandidate> MuonAnalysisb;
DEFINE_FWK_MODULE(MuonAnalysisb);

