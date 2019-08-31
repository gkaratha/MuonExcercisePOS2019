import FWCore.ParameterSet.Config as cms

muonGenMatching =  cms.EDAnalyzer('MuonGenMatching',
    muonSrc = cms.InputTag('slimmedMuons'),
    vertexSrc = cms.InputTag('selectedPrimaryVertices'),
    genParticleSrc = cms.InputTag('prunedGenParticles'),
)
