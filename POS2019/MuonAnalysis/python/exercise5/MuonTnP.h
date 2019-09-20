//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 11 18:01:22 2019 by ROOT version 6.16/00
// from TTree mytree/mytree
// found on file: DYJetsTree_1.root
//////////////////////////////////////////////////////////

#ifndef MuonTnP_h
#define MuonTnP_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "iostream"
using std::vector;

class MuonTnP {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           run_number;
   Int_t           ls;
   Float_t         pvertex_x;
   Float_t         pvertex_y;
   Float_t         pvertex_z;
   Float_t         bspot_x;
   Float_t         bspot_y;
   Float_t         bspot_z;
   vector<float>   *genmuon_pt;
   vector<float>   *genmuon_eta;
   vector<float>   *genmuon_phi;
   vector<float>   *genmuon_charge;
   vector<float>   *genmuon_motmId;
   vector<float>   *genmuon_grandmomId;
   Bool_t          hltpath;
   vector<float>   *l1muon_pt;
   vector<float>   *l1muon_eta;
   vector<float>   *l1muon_phi;
   vector<float>   *hltmuon_pt;
   vector<float>   *HLTmuon_eta;
   vector<float>   *hltmuon_phi;
   Int_t           nmuon;
   vector<float>   *muon_pt;
   vector<float>   *muon_eta;
   vector<float>   *muon_phi;
   vector<float>   *muon_charge;
   vector<float>   *muon_vx;
   vector<float>   *muon_vy;
   vector<float>   *muon_vz;
   vector<float>   *muon_iso;
   vector<bool>    *muon_soft;
   vector<bool>    *muon_loose;
   vector<bool>    *muon_medium;
   vector<bool>    *muon_tight;
   Int_t           njpsi;
   vector<float>   *jpsi_pt;
   vector<float>   *jpsi_eta;
   vector<float>   *jpsi_phi;
   vector<float>   *jpsi_ctxy;
   vector<float>   *jpsi_x;
   vector<float>   *jpsi_y;
   vector<float>   *jpsi_z;
   vector<float>   *jpsi_prob;
   vector<float>   *jpsi_mass;
   vector<float>   *jpsi_mass_unfit;
   vector<float>   *jpsi_mu1pt;
   vector<float>   *jpsi_mu1eta;
   vector<float>   *jpsi_mu1phi;
   vector<float>   *jpsi_mu2pt;
   vector<float>   *jpsi_mu2eta;
   vector<float>   *jpsi_mu2phi;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run_number;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_pvertex_x;   //!
   TBranch        *b_pvertex_y;   //!
   TBranch        *b_pvertex_z;   //!
   TBranch        *b_bspot_x;   //!
   TBranch        *b_bspot_y;   //!
   TBranch        *b_bspot_z;   //!
   TBranch        *b_genmuon_pt;   //!
   TBranch        *b_genmuon_eta;   //!
   TBranch        *b_genmuon_phi;   //!
   TBranch        *b_genmuon_charge;   //!
   TBranch        *b_genmuon_motmId;   //!
   TBranch        *b_genmuon_grandmomId;   //!
   TBranch        *b_hltpath;   //!
   TBranch        *b_l1muon_pt;   //!
   TBranch        *b_l1muon_eta;   //!
   TBranch        *b_l1muon_phi;   //!
   TBranch        *b_hltmuon_pt;   //!
   TBranch        *b_HLTmuon_eta;   //!
   TBranch        *b_hltmuon_phi;   //!
   TBranch        *b_nmuon;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_vx;   //!
   TBranch        *b_muon_vy;   //!
   TBranch        *b_muon_vz;   //!
   TBranch        *b_muon_iso;   //!
   TBranch        *b_muon_soft;   //!
   TBranch        *b_muon_loose;   //!
   TBranch        *b_muon_medium;   //!
   TBranch        *b_muon_tight;   //!
   TBranch        *b_njpsi;   //!
   TBranch        *b_jpsi_pt;   //!
   TBranch        *b_jpsi_eta;   //!
   TBranch        *b_jpsi_phi;   //!
   TBranch        *b_jpsi_ctxy;   //!
   TBranch        *b_jpsi_x;   //!
   TBranch        *b_jpsi_y;   //!
   TBranch        *b_jpsi_z;   //!
   TBranch        *b_jpsi_prob;   //!
   TBranch        *b_jpsi_mass;   //!
   TBranch        *b_jpsi_mass_unfit;   //!
   TBranch        *b_jpsi_mu1pt;   //!
   TBranch        *b_jpsi_mu1eta;   //!
   TBranch        *b_jpsi_mu1phi;   //!
   TBranch        *b_jpsi_mu2pt;   //!
   TBranch        *b_jpsi_mu2eta;   //!
   TBranch        *b_jpsi_mu2phi;   //!

   MuonTnP(TTree *tree=0);
   virtual ~MuonTnP();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(float tagPt, float probePt, float etaCut, float mLow, float mUp, float tagIso, float probeIso, int maxEvts, TString outFile, bool isJpsi = false);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


};

#endif

#ifdef MuonTnP_cxx
MuonTnP::MuonTnP(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/SingleMuon/SingleMuon_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/SingleMuon/SingleMuon_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/SingleMuon/SingleMuon_1.root:/demo");
      dir->GetObject("mytree",tree);

   }
   Init(tree);
}

MuonTnP::~MuonTnP()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MuonTnP::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MuonTnP::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MuonTnP::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genmuon_pt = 0;
   genmuon_eta = 0;
   genmuon_phi = 0;
   genmuon_charge = 0;
   genmuon_motmId = 0;
   genmuon_grandmomId = 0;
   l1muon_pt = 0;
   l1muon_eta = 0;
   l1muon_phi = 0;
   hltmuon_pt = 0;
   HLTmuon_eta = 0;
   hltmuon_phi = 0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_charge = 0;
   muon_vx = 0;
   muon_vy = 0;
   muon_vz = 0;
   muon_iso = 0;
   muon_soft = 0;
   muon_loose = 0;
   muon_medium = 0;
   muon_tight = 0;
   jpsi_pt = 0;
   jpsi_eta = 0;
   jpsi_phi = 0;
   jpsi_ctxy = 0;
   jpsi_x = 0;
   jpsi_y = 0;
   jpsi_z = 0;
   jpsi_prob = 0;
   jpsi_mass = 0;
   jpsi_mass_unfit = 0;
   jpsi_mu1pt = 0;
   jpsi_mu1eta = 0;
   jpsi_mu1phi = 0;
   jpsi_mu2pt = 0;
   jpsi_mu2eta = 0;
   jpsi_mu2phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run_number", &run_number, &b_run_number);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("pvertex_x", &pvertex_x, &b_pvertex_x);
   fChain->SetBranchAddress("pvertex_y", &pvertex_y, &b_pvertex_y);
   fChain->SetBranchAddress("pvertex_z", &pvertex_z, &b_pvertex_z);
   fChain->SetBranchAddress("bspot_x", &bspot_x, &b_bspot_x);
   fChain->SetBranchAddress("bspot_y", &bspot_y, &b_bspot_y);
   fChain->SetBranchAddress("bspot_z", &bspot_z, &b_bspot_z);
   fChain->SetBranchAddress("genmuon_pt", &genmuon_pt, &b_genmuon_pt);
   fChain->SetBranchAddress("genmuon_eta", &genmuon_eta, &b_genmuon_eta);
   fChain->SetBranchAddress("genmuon_phi", &genmuon_phi, &b_genmuon_phi);
   fChain->SetBranchAddress("genmuon_charge", &genmuon_charge, &b_genmuon_charge);
   fChain->SetBranchAddress("genmuon_motmId", &genmuon_motmId, &b_genmuon_motmId);
   fChain->SetBranchAddress("genmuon_grandmomId", &genmuon_grandmomId, &b_genmuon_grandmomId);
   fChain->SetBranchAddress("hltpath", &hltpath, &b_hltpath);
   fChain->SetBranchAddress("l1muon_pt", &l1muon_pt, &b_l1muon_pt);
   fChain->SetBranchAddress("l1muon_eta", &l1muon_eta, &b_l1muon_eta);
   fChain->SetBranchAddress("l1muon_phi", &l1muon_phi, &b_l1muon_phi);
   fChain->SetBranchAddress("hltmuon_pt", &hltmuon_pt, &b_hltmuon_pt);
   fChain->SetBranchAddress("HLTmuon_eta", &HLTmuon_eta, &b_HLTmuon_eta);
   fChain->SetBranchAddress("hltmuon_phi", &hltmuon_phi, &b_hltmuon_phi);
   fChain->SetBranchAddress("nmuon", &nmuon, &b_nmuon);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_charge", &muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_vx", &muon_vx, &b_muon_vx);
   fChain->SetBranchAddress("muon_vy", &muon_vy, &b_muon_vy);
   fChain->SetBranchAddress("muon_vz", &muon_vz, &b_muon_vz);
   fChain->SetBranchAddress("muon_iso", &muon_iso, &b_muon_iso);
   fChain->SetBranchAddress("muon_soft", &muon_soft, &b_muon_soft);
   fChain->SetBranchAddress("muon_loose", &muon_loose, &b_muon_loose);
   fChain->SetBranchAddress("muon_medium", &muon_medium, &b_muon_medium);
   fChain->SetBranchAddress("muon_tight", &muon_tight, &b_muon_tight);
   fChain->SetBranchAddress("njpsi", &njpsi, &b_njpsi);
   fChain->SetBranchAddress("jpsi_pt", &jpsi_pt, &b_jpsi_pt);
   fChain->SetBranchAddress("jpsi_eta", &jpsi_eta, &b_jpsi_eta);
   fChain->SetBranchAddress("jpsi_phi", &jpsi_phi, &b_jpsi_phi);
   fChain->SetBranchAddress("jpsi_ctxy", &jpsi_ctxy, &b_jpsi_ctxy);
   fChain->SetBranchAddress("jpsi_x", &jpsi_x, &b_jpsi_x);
   fChain->SetBranchAddress("jpsi_y", &jpsi_y, &b_jpsi_y);
   fChain->SetBranchAddress("jpsi_z", &jpsi_z, &b_jpsi_z);
   fChain->SetBranchAddress("jpsi_prob", &jpsi_prob, &b_jpsi_prob);
   fChain->SetBranchAddress("jpsi_mass", &jpsi_mass, &b_jpsi_mass);
   fChain->SetBranchAddress("jpsi_mass_unfit", &jpsi_mass_unfit, &b_jpsi_mass_unfit);
   fChain->SetBranchAddress("jpsi_mu1pt", &jpsi_mu1pt, &b_jpsi_mu1pt);
   fChain->SetBranchAddress("jpsi_mu1eta", &jpsi_mu1eta, &b_jpsi_mu1eta);
   fChain->SetBranchAddress("jpsi_mu1phi", &jpsi_mu1phi, &b_jpsi_mu1phi);
   fChain->SetBranchAddress("jpsi_mu2pt", &jpsi_mu2pt, &b_jpsi_mu2pt);
   fChain->SetBranchAddress("jpsi_mu2eta", &jpsi_mu2eta, &b_jpsi_mu2eta);
   fChain->SetBranchAddress("jpsi_mu2phi", &jpsi_mu2phi, &b_jpsi_mu2phi);

   std::cout << "Loaded TTree with " << fChain->GetEntries() << " entries\n";
   Notify();
}

Bool_t MuonTnP::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MuonTnP::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MuonTnP::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MuonTnP_cxx
