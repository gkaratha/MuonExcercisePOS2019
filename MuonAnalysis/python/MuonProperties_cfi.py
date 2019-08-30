import FWCore.ParameterSet.Config as cms

muonProperties =  cms.EDAnalyzer('MuonProperties',
    muonSrc = cms.InputTag('slimmedMuons'),
    vertexSrc = cms.InputTag('offlineSlimmedPrimaryVertices'),
    genParticleSrc = cms.InputTag('prunedGenParticles'),
    isQCD = cms.untracked.bool(False)
)
