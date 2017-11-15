import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outFilename', 'MiniEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('inputFormat', 'PAT',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "format of the input files (PAT or RECO)"
                 )
options.register('skim', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "skim events with one lepton and 2 jets"
                 )
options.register('updateJEC', '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "Name of the SQLite file (with path and extension) used to update the jet collection to the latest JEC and the era of the new JEC"
                )
options.parseArguments()

process = cms.Process("MiniAnalysis")

# Geometry, GT, and other standard sequences
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '91X_upgrade2023_realistic_v3', '')

# Log settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('MyAna')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
)

# Input
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3) ) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17MiniAOD/DiPhotonJetsBox_MGG-80toInf_14TeV-Sherpa/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v1/150000/E01EFC0F-13B7-E711-A90B-FA163E4C681C.root'
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17MiniAOD/VBFHToGG_M125_14TeV_amcatnlo_pythia8/MINIAODSIM/noPU_93X_upgrade2023_realistic_v2-v1/150000/40333064-FDBB-E711-ADFD-0025905C2CD2.root'
        'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17MiniAOD/VBFHToGG_M125_14TeV_amcatnlo_pythia8/MINIAODSIM/noPU_93X_upgrade2023_realistic_v2-v1/150000/661130DC-21BC-E711-939D-3417EBE2EE2D.root'
    ),
    secondaryFileNames = cms.untracked.vstring(
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DiPhotonJetsBox_MGG-80toInf_14TeV-Sherpa/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/24BB7BB2-A9B4-E711-9DC9-FA163E7FFB3C.root')
       'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/VBFHToGG_M125_14TeV_amcatnlo_pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/150000/00C13FEA-40B3-E711-857C-008CFAFBE5E0.root','root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/VBFHToGG_M125_14TeV_amcatnlo_pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/150000/D2CB12B4-3FB3-E711-A6DE-3417EBE2F319.root','root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/VBFHToGG_M125_14TeV_amcatnlo_pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/150000/E0B60B0C-A5B3-E711-AC36-008CFAF71768.root','root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/VBFHToGG_M125_14TeV_amcatnlo_pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/150000/0CA2BC2F-A5B3-E711-9E2B-008CFAFBF576.root','root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/VBFHToGG_M125_14TeV_amcatnlo_pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/150000/148522BA-67BB-E711-9D2B-3417EBE64567.root') 
)
if (options.inputFormat.lower() == "reco"):
    process.source.fileNames = cms.untracked.vstring(*(
        '/store/mc/PhaseIITDRSpring17DR/TTToSemiLepton_TuneCUETP8M1_14TeV-powheg-pythia8/AODSIM/PU200_91X_upgrade2023_realistic_v3-v1/120000/000CD008-7A58-E711-82DB-1CB72C0A3A61.root',
    ))

# HGCAL Photon ID
process.load("RecoEgamma.Phase2InterimID.phase2EgammaPAT_cff")
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True)
)

# Get new JEC from an SQLite file rather than a GT
if options.updateJEC:
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string('sqlite_file:'+options.updateJEC[0]),
                               toGet =  cms.VPSet(
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                     tag = cms.string("JetCorrectorParametersCollection_"+options.updateJEC[1]+"_AK4PFPuppi"),
                     label = cms.untracked.string("AK4PFPuppi"))
            )
                               )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")

process.source.inputCommands = cms.untracked.vstring("keep *")

# Pre-skim weight counter
process.weightCounter = cms.EDAnalyzer('WeightCounter')

# Skim filter
muonLabel = "slimmedMuons"
elecLabel = "slimmedElectrons"
photLabel = "phase2Photons"
if (options.inputFormat.lower() == "reco"):
    photLabel = "slimmedPhotons"
if options.updateJEC:
    jetLabel = "updatedPatJetsUpdatedJECAK4PFPuppi"
else:    
    jetLabel = "slimmedJets"
if (options.inputFormat.lower() == "reco"):
    muonLabel = "muons"
    elecLabel = "ecalDrivenGsfElectrons"
    if options.updateJEC:
        jetLabel = "ak4PUPPIJetsL1FastL2L3"
    else:    
        jetLabel = "ak4PUPPIJets"
process.selectedMuons = cms.EDFilter("CandPtrSelector",
                                     src = cms.InputTag(muonLabel),
                                     cut = cms.string("pt>10 && abs(eta)<3")
                                     )
process.selectedElectrons = cms.EDFilter("CandPtrSelector",
                                         src = cms.InputTag(elecLabel),
                                         cut = cms.string("pt>10 && abs(eta)<3")
                                         )
process.selectedPhotons = cms.EDFilter("CandPtrSelector",
                                         src = cms.InputTag(photLabel),
                                         cut = cms.string("pt>10 && abs(eta)<3")
                                         )
process.selectedJets = cms.EDFilter("CandPtrSelector",
                                    src = cms.InputTag(jetLabel),
                                    cut = cms.string("pt>20 && abs(eta)<5")
                                    )
process.allLeps = cms.EDProducer("CandViewMerger",
                                 src = cms.VInputTag(
                                                     cms.InputTag("selectedElectrons"),
                                                     cms.InputTag("selectedMuons")
                                                     )
                                 )
process.countLeps = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("allLeps"),
                                 minNumber = cms.uint32(1)
                                 )
process.countPhotons = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("selectedPhotons"),
                                 minNumber = cms.uint32(0)
                                 )
process.countJets = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("selectedJets"),
                                 minNumber = cms.uint32(2)
                                 )
process.preYieldFilter = cms.Sequence(process.selectedMuons+process.selectedElectrons+process.allLeps+process.countLeps+process.selectedPhotons+process.countPhotons+process.selectedJets+process.countJets)


# run Puppi 
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.load('CommonTools/PileupAlgos/PhotonPuppi_cff')
process.load('CommonTools/PileupAlgos/softKiller_cfi')
from CommonTools.PileupAlgos.PhotonPuppi_cff        import setupPuppiPhoton
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppies
makePuppies(process)
process.particleFlowNoLep = cms.EDFilter("PdgIdCandViewSelector",
                                    src = cms.InputTag("particleFlow"), 
                                    pdgId = cms.vint32( 1,2,22,111,130,310,2112,211,-211,321,-321,999211,2212,-2212 )
                                    )
process.puppiNoLep = process.puppi.clone(candName = cms.InputTag('particleFlowNoLep'))
process.load("PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi")
process.load("PhysicsTools.PatAlgos.slimming.offlineSlimmedPrimaryVertices_cfi")
process.load("PhysicsTools.PatAlgos.slimming.packedPFCandidates_cfi")

# recluster jets
process.load('RecoJets/Configuration/RecoPFJets_cff')
process.ak4PUPPIJets  = process.ak4PFJets.clone(rParam=0.4, src = cms.InputTag('puppi'))

# recompute MET
process.load('RecoMET.METProducers.PFMET_cfi')
process.puppiMet = process.pfMet.clone()
process.puppiMet.src = cms.InputTag('puppi')

# PF cluster producer for HFCal ID
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cff")

# jurassic track isolation
# https://indico.cern.ch/event/27568/contributions/1618615/attachments/499629/690192/080421.Isolation.Update.RecHits.pdf
process.load("RecoEgamma.EgammaIsolationAlgos.electronTrackIsolationLcone_cfi")
process.electronTrackIsolationLcone.electronProducer = cms.InputTag("ecalDrivenGsfElectrons")
process.electronTrackIsolationLcone.intRadiusBarrel = 0.04
process.electronTrackIsolationLcone.intRadiusEndcap = 0.04

# analysis
moduleName = "MiniFromPat"    
if (options.inputFormat.lower() == "reco"):
    moduleName = "MiniFromReco"
process.ntuple = cms.EDAnalyzer(moduleName)
process.load("PhaseTwoAnalysis.NTupler."+moduleName+"_cfi")
if (options.inputFormat.lower() == "reco"):
    process.ntuple.pfCandsNoLep = "puppiNoLep"
    process.ntuple.met = "puppiMet"
    if options.updateJEC:
        # This will load several ESProducers and EDProducers which make the corrected jet collections
        # In this case the collection will be called ak4PUPPIJetsL1FastL2L3
        process.load('PhaseTwoAnalysis.Jets.JetCorrection_cff')
        process.ntuple.jets = "ak4PUPPIJetsL1FastL2L3"
    else:
        # This simply switches the default AK4PFJetsCHS collection to the ak4PUPPIJets collection now that it has been produced
        process.ntuple.jets = "ak4PUPPIJets"
else:
    if options.updateJEC:
        # The updateJetCollection function will uncorred the jets from MiniAOD and then recorrect them using the current
        #  set of JEC in the event setup
        # The new name of the updated jet collection becomes updatedPatJetsUpdatedJECAK4PFPuppi
        from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
        updateJetCollection(process,
                            jetSource = cms.InputTag('slimmedJetsPuppi'),
                            postfix = 'UpdatedJECAK4PFPuppi',
                            jetCorrections = ('AK4PFPuppi', ['L1FastJet','L2Relative','L3Absolute'], 'None')
                            )
        process.ntuple.jets = "updatedPatJetsUpdatedJECAK4PFPuppi"

# output
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                   )

# run
if (options.inputFormat.lower() == "reco"):
    process.puSequence = cms.Sequence(process.primaryVertexAssociation * process.pfNoLepPUPPI * process.puppi * process.particleFlowNoLep * process.puppiNoLep * process.offlineSlimmedPrimaryVertices * process.packedPFCandidates * process.muonIsolationPUPPI * process.muonIsolationPUPPINoLep * process.ak4PUPPIJets * process.puppiMet)

if options.skim:
    if (options.inputFormat.lower() == "reco"):
        if options.updateJEC:
            process.p = cms.Path(process.weightCounter * process.phase2Egamma * process.electronTrackIsolationLcone * process.particleFlowRecHitHGCSeq * process.puSequence * process.ak4PFPuppiL1FastL2L3CorrectorChain * process.ak4PUPPIJetsL1FastL2L3 * process.preYieldFilter * process.ntuple)
        else:
            process.p = cms.Path(process.weightCounter * process.phase2Egamma * process.electronTrackIsolationLcone * process.particleFlowRecHitHGCSeq * process.puSequence * process.preYieldFilter * process.ntuple)
    else:
        if options.updateJEC:
            process.p = cms.Path(process.weightCounter*process.phase2Egamma*process.preYieldFilter*process.patJetCorrFactorsUpdatedJECAK4PFPuppi * process.updatedPatJetsUpdatedJECAK4PFPuppi * process.ntuple)
        else:
            process.p = cms.Path(process.weightCounter*process.phase2Egamma*process.preYieldFilter*process.ntuple)
else:
    if (options.inputFormat.lower() == "reco"):
        if options.updateJEC:
            process.p = cms.Path(process.weightCounter*process.electronTrackIsolationLcone * process.particleFlowRecHitHGCSeq * process.puSequence * process.ak4PFPuppiL1FastL2L3CorrectorChain * process.ak4PUPPIJetsL1FastL2L3 * process.phase2Egamma * process.ntuple)
        else:
            process.p = cms.Path(process.weightCounter*process.electronTrackIsolationLcone * process.particleFlowRecHitHGCSeq * process.puSequence * process.phase2Egamma * process.ntuple)
    else:
        if options.updateJEC:
            process.p = cms.Path(process.weightCounter*process.patJetCorrFactorsUpdatedJECAK4PFPuppi * process.updatedPatJetsUpdatedJECAK4PFPuppi * process.phase2Egamma * process.ntuple)
	else:    
            process.p = cms.Path(process.weightCounter*process.phase2Egamma*process.ntuple)
