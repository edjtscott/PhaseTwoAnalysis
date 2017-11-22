from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = 'VBFHToGG_M125_14TeV_amcatnlo_pythia8-PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1'
#config.General.requestName = 'VBFHToGG_M125_14TeV_amcatnlo_pythia8--PU200_93X_upgrade2023_realistic_v2-v2_2ndAttempt'
#config.General.requestName = 'ggHToGG_M125_14TeV_amcatnloFXFX_pythia8--noPU_93X_upgrade2023_realistic_v2-v1_2ndAttempt'
config.General.requestName = 'ggHToGG_M125_14TeV_amcatnloFXFX_pythia8--PU200_93X_upgrade2023_realistic_v2-v2_2ndAttempt'
config.General.workArea = 'crab_tasks/'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'scripts/produceNtuples_cfg.py'
config.JobType.pyCfgParams= ['skim=False','inputFormat=PAT','outFilename=MiniEvents.root']
# Uncomment the following line when running on PAT events
config.JobType.inputFiles = ['PhaseIIFall17_V3_MC.db']
config.JobType.outputFiles = ['MiniEvents.root']

config.section_("Data")
#config.Data.inputDataset = '/VBFHToGG_M125_14TeV_amcatnlo_pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'
#config.Data.inputDataset = '/VBFHToGG_M125_14TeV_amcatnlo_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
#config.Data.inputDataset = '/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/MINIAODSIM'
config.Data.inputDataset = '/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/PhaseIITDRFall17MiniAOD-PU200_93X_upgrade2023_realistic_v2-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 5
config.Data.unitsPerJob = 1
# Uncomment to run on a fraction of the dataset
#config.Data.totalUnits = 5
#config.Data.outLFNDirBase = '/store/group/dpg_hgcal/comm_hgcal/escott/VBF_PU0_17Nov17/'
#config.Data.outLFNDirBase = '/store/group/dpg_hgcal/comm_hgcal/escott/VBF_PU200_20Nov17/'
#config.Data.outLFNDirBase = '/store/group/dpg_hgcal/comm_hgcal/escott/ggH_PU0_20Nov17/'
config.Data.outLFNDirBase = '/store/group/dpg_hgcal/comm_hgcal/escott/ggH_PU200_20Nov17/'
config.Data.publication = False

config.Data.useParent = True # need to run on GEN-SIM-RECO to apply photon ID

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.ignoreGlobalBlacklist = True
