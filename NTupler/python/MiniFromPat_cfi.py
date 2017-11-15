import FWCore.ParameterSet.Config as cms

ntuple = cms.EDAnalyzer('MiniFromPat',
        pileup        = cms.uint32(200),
        vertices      = cms.InputTag("offlineSlimmedPrimaryVertices"),
        electrons     = cms.InputTag("slimmedElectrons"),
        photons       = cms.InputTag("phase2Photons"),
        beamspot      = cms.InputTag("offlineBeamSpot"),
        conversions   = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
        muons         = cms.InputTag("slimmedMuons"),
        jets          = cms.InputTag("slimmedJetsPuppi"),
        mets          = cms.InputTag("slimmedMETsPuppi"),
        genParts      = cms.InputTag("packedGenParticles"),
        genVertices   = cms.InputTag("genParticles", "xyz0", "HLT"),
        genJets       = cms.InputTag("slimmedGenJets"),
)
