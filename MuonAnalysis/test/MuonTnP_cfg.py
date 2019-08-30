import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('TnPAnalysis')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/user/r/roy4/data/E6284925-C335-E711-B792-1CB72C1B2EF4.root'),
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

# Primary Vertex Selector
process.selectedPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlineSlimmedPrimaryVertices'),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && abs(position.Rho) <= 2"),
  filter = cms.bool(True)                                          
)


process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('muonTnP.root')
)

process.load('MuonPOGLongExercise.MuonTnP.MuonTnP_cfi')
                    
process.p = cms.Path(process.selectedPrimaryVertices*process.muonTnP)

