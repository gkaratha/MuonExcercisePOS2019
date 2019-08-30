import FWCore.ParameterSet.Config as cms

muonTnP =  cms.EDAnalyzer('MuonTnP',
    muonSrc = cms.InputTag('slimmedMuons'),
    vertexSrc = cms.InputTag('selectedPrimaryVertices'),
    hltTag = cms.InputTag('TriggerResults','','HLT'),
    triggerObjectTag = cms.InputTag('selectedPatTrigger'),
)
