from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters

idx = '00000'
CFG = 'QCD_Pt_80_170_v2_%s'%idx

# To submit to crab:
# crab submit -c crabConfig_data.py
# To check job status:
# crab status -d <config.General.workArea>/<config.General.requestName># To resubmit jobs:
# crab resubmit -d <config.General.workArea>/<config.General.requestName>

# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>
config.General.workArea = 'crab_MC'
config.General.requestName = CFG
config.General.transferOutputs = True
config.General.transferLogs = False

# CMS cfg file goes here:
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'MLAnalyzer/RecHitAnalyzer/python/ConfFile_cfg.py' # analyzer cfg file
#config.JobType.maxMemoryMB = 2800

# Define input and units per job here:
config.Data.userInputFiles = open('MLAnalyzer/LISTS/CMS_MonteCarlo2012_Summer12_DR53X_QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6_AODSIM_PU_RD1_START53_V7N-v3_%s_file_index.txt'%idx).readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30 # units: as defined by config.Data.splitting
config.Data.totalUnits = -1 # -1: all inputs. total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
#config.Data.totalUnits = 10 # test production
config.Data.publication = False

# Output files will be stored in config.Site.storageSite at directory:
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/lpcml/' # add your username as subdirectory
config.Data.outputPrimaryDataset = 'Run2017_v2'
config.Data.outputDatasetTag = config.General.requestName
