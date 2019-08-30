import FWCore.ParameterSet.Config as cms

muonTnPmctruth =  cms.EDAnalyzer('MuonTnPMCtruth',
    muonSrc = cms.InputTag('slimmedMuons'),
    vertexSrc = cms.InputTag('selectedPrimaryVertices'),
    genParticleSrc = cms.InputTag('prunedGenParticles')
)
