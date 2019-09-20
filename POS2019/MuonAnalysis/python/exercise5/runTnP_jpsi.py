#!/usr/bin/env python
import sys
from ROOT import *


gROOT.ProcessLine(".L MuonTnP.C++")
chain= TChain("demo/mytree","")
outFile=""
if sys.argv[1] ==  "MC":
    #dpsiMC
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BtoJpsiK_MC/output_1.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BtoJpsiK_MC/output_2.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BtoJpsiK_MC/output_3.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BtoJpsiK_MC/output_4.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BtoJpsiK_MC/output_5.root")    
    outFile = "tnp_jpsi_mc.root"
else:
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_1.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_10.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_11.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_12.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_13.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_14.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_15.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_16.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_17.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_18.root")
    chain.Add("/eos/user/r/roy4/data/cmspos2019/muPOG/trees/BParking_Run2018B_M_2p5_4/output_19.root")
    outFile = "tnp_jpsi_data.root"


tnp=MuonTnP(chain)
# Loop(tagPt, probePt, etaCut, mLow, mUp, tagIso,  probeIso, int maxEvts, TString outFile, bool isJpsi);
#For Z selection
tnp.Loop(10., 2., 2., 3., 3.15, 9999., 0.2, 2000000,  outFile, 1)
