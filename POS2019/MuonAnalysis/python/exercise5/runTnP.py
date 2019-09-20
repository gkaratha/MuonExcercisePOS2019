#!/usr/bin/env python
import sys
from ROOT import *


gROOT.ProcessLine(".L MuonTnP.C++")
chain= TChain("demo/mytree","")
outFile=""
if sys.argv[1] ==  "MC":
    #dyjetMC
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/dyjets/dyjetsMerged.root")
    outFile = "tnp_z_mc.root"
else:
    ##SingleMuondata
    #chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/SingleMuon/SingleMuon_1.root")
    #chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/SingleMuon/SingleMuon_2.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/SingleMuon/SingleMuon_3.root")
    outFile = "tnp_z_data.root"


tnp=MuonTnP(chain)
# Loop(tagPt, probePt, etaCut, mLow, float mUp, float tagIso, float probeIso, int maxEvts, TString outFile);
#For Z selection
tnp.Loop(26., 20., 2.4, 85., 95., 0.2, 0.2, 2000000,  outFile)
