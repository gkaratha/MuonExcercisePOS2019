#!/usr/bin/env python
from ROOT import *
import sys

def makeEffhistosfromlist(inFile, inputhlist, effplots, var):
    for idx in range(1, len(inputhlist)):
        print 'Taking ratio of',  inputhlist[idx], inputhlist[idx-1]
        hden=tnpHfile.Get(inputhlist[idx-1])
        hnum=tnpHfile.Get(inputhlist[idx])
        gae_pt=TGraphAsymmErrors()
        gae_pt.SetName('eff_' + inputhlist[idx])
        gae_pt.Divide(hnum, hden, "cl=0.683 b(1,1) mode")
        gae_pt.GetXaxis().SetTitle(var)
        effplots.append(gae_pt)
 


tnpHfile=TFile.Open(sys.argv[1])

histosPt=['hmuPt_probe_den', 'hmuPt_probe_num_looseId', 'hmuPt_probe_num_mediumId', 'hmuPt_probe_num_tightId', 'hmuPt_probe_num_iso']
histosEta=['hmuEta_probe_den', 'hmuEta_probe_num_looseId', 'hmuEta_probe_num_mediumId', 'hmuEta_probe_num_tightId', 'hmuEta_probe_num_iso']

effplots=[]

makeEffhistosfromlist(tnpHfile, histosPt, effplots, 'p_{T} GeV') 
makeEffhistosfromlist(tnpHfile, histosEta, effplots, '#eta') 

outFile=TFile.Open('effPlots.root', 'recreate')
outFile.cd()
for h in effplots:
    h.Write()
outFile.Save()
outFile.Close()
tnpHfile.Close()
