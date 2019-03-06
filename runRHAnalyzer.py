import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/eusai/step2_ttbar_p8_03/step2_qcd8_1109.root'
#inputFiles_='file:../step_full.root'
#inputFiles_='file:/uscms/home/bburkle/nobackup/working_area/CMSSW_5_3_32/src/opendatadnn/step2_test.root'
#inputFiles_='file:/uscms/home/bburkle/nobackup/working_area/CMSSW_5_3_32/src/MLAnalyzer/test/step2_ttbarOD.root'
#inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/eusai/CRAB_UserFiles/step2_QCD600to3000_01/190213_183439/0000/step2_QCDPt_15_3000_Flat_V27_961.root'
inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/eusai/CRAB_UserFiles/step2_TtbarFromOpen_01/190220_154334/0000/step2_TtbarFromOpen_101.root'

isTTbar_ = 1

#maxEvents_=100
#skipEvents_=0#
#outputFile_ = 'test.root'
outputFile_ = 'test/ttbar_new-production_test.root'

cmd="cmsRun %s inputFiles=%s outputFile=%s isTTbar=%d" %(cfg,inputFiles_,outputFile_,isTTbar_)
print '%s'%cmd
os.system(cmd)
