import os
#from glob import glob
import re

#listname = 'CMS_MonteCarlo2012_Summer12_DR53X_DiPhotonBorn_Pt-25To250_8TeV_ext-pythia6_AODSIM_PU_RD1_START53_V7N-v2_10000_file_index.txt'
#listname = 'CMS_MonteCarlo2012_Summer12_DR53X_DiPhotonBox_Pt-25To250_8TeV_ext-pythia6_AODSIM_PU_RD1_START53_V7N-v2_10000_file_index.txt'
#listname = 'CMS_MonteCarlo2012_Summer12_DR53X_GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV_ext-pythia6_AODSIM_PU_RD1_START53_V7N-v1_20000_file_index.txt'
#listname = 'CMS_MonteCarlo2012_Summer12_DR53X_GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV_ext-pythia6_AODSIM_PU_RD1_START53_V7N-v1_20001_file_index.txt'
#listname = 'CMS_MonteCarlo2012_Summer12_DR53X_GluGluHToGG_M-125_8TeV-pythia6_AODSIM_PU_RD1_START53_V7N-v2_00000_file_index.txt'
#listname = 'CMS_MonteCarlo2012_Summer12_DR53X_GluGluHToGG_M-125_8TeV-pythia6_AODSIM_PU_RD1_START53_V7N-v2_10000_file_index.txt'
#listname = 'CMS_MonteCarlo2012_Summer12_DR53X_GluGluHToGG_M-125_8TeV-pythia6_AODSIM_PU_RD1_START53_V7N-v2_20000_file_index.txt'
#listname = 'CMS_MonteCarlo2012_Summer12_DR53X_WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_AODSIM_PU_RD1_START53_V7N-v1_10000_file_index.txt'
#listname = 'CMS_MonteCarlo2012_Summer12_DR53X_QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6_AODSIM_PU_RD1_START53_V7N-v3_00000_file_index.txt'
listname = 'topgun_mtop85To250GeV_pT10To800GeV_eta1p4_noPU_noTkrh_file_list.txt'

#eosDir='/eos/cms/store/user/mandrews/OPENDATA'
#decay='%s_FEVTDEBUG'%args.decay
#

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_ = ['file:%s'%path for path in glob('%s/FEVTDEBUG/%s/*/*/step*root'%(eosDir,decay))]
#
#listname = 'list_%s.txt'%decay
#with open(listname, 'w') as list_file:
#    for inputFile in inputFiles_:
#        list_file.write("%s\n" % inputFile)
#listname = '%s.txt'%args.decay
#proc = re.search('CMS_MonteCarlo2012_Summer12_DR53X_(.+?)_', listname).group(1)
#idx = re.search('-v.+?_(.+?)_file_index', listname).group(1)
#decay = '%s_%s'%(proc,idx)
#print(decay)

maxEvents_=-1
maxEvents_=100000
skipEvents_=0

isTTbar_ = 1
outputFile_='/afs/cern.ch/user/m/mandrews/public/ForDarya/topgun_ntuple.root'

#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
#cmd="cmsRun %s inputFiles_load=LISTS/%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,listname,maxEvents_,skipEvents_,outfile)
cmd="cmsRun %s inputFiles_load=LISTS/%s outputFile=%s isTTbar=%d maxEv=%d" %(cfg,listname,outputFile_,isTTbar_,maxEvents_)
#print '%s'%cmd
os.system(cmd)

#os.system('mv cEB*.eps %s/'%(inputTag))
