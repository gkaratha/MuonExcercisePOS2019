#!/usr/bin/env python
from ROOT import *


gROOT.ProcessLine(".L MuonTnP.C++")
chain= TChain("demo/mytree","")
chain.Add("dyjets_1.root")
chain.Add("dyjets_2.root")
chain.Add("dyjets_3.root")

tnp=MuonTnP()
tnp.Loop(24., 20., 10., 0.15, 0.15)
