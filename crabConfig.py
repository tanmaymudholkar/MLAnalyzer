from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#dataset_ = '/RSGluonToTT_M-4000_TuneCUETP8M1_14TeV-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
dataset_ = '/TT_Mtt1500toInf_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
#dataset_ = '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
outputFile_ = dataset_.split('/')[1]+'.root'
maxEvents_ = -1
skipEvents_ = 0


#the_name = 'RS_Gluon_Image_Ntuplizer'
the_name = 'TTbar_Image_NTuplizer_1f'
#the_name = 'QCD_Background_Ntuplizer'
config.General.requestName = the_name
config.General.workArea = '/uscms/home/bburkle/nobackup/working_area/CMSSW_9_3_0/src/MLAnalyzer'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RecHitAnalyzer/python/ConfFile_cfg.py'
config.JobType.pyCfgParams = ['outputFile=%s'%outputFile_]
config.JobType.maxMemoryMB = 4000
#config.JobType.pyCfgParams = ['maxEvents=%s'%maxEvents_,'outputFile=%s'%outputFile_,'skipEvents=%s'%skipEvents_]
#config.JobType.pluginName = 'RecHitAnalyzer/python/ConfFile_cfg.py outputFile=%s maxEvents=%s skipEvents=%s' % (outputFile_, maxEvents_, skipEvents_)

config.Data.inputDataset = dataset_
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 #if splitting by files
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/E2E/root_files/' % (getUsernameFromSiteDB())
config.Data.outputDatasetTag = the_name

config.Site.storageSite = 'T3_US_FNALLPC'
