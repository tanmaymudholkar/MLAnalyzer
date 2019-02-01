import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/FEVTDEBUG/h24gamma_1j_10K_100MeV_FEVTDEBUG_2016_25ns_Moriond17MC_PoissonOOTPU/180109_112606/0000/step_full_1.root'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/AODSIM/T7WgStealth_800_200_AODSIM_noPU/180227_224514/0000/step_AODSIM_1.root'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/AODSIM/WZToJets_TuneCUETP8M1_13TeV_pythia8_noPU_AODSIM/0/0000/step_full_1.root'
#inputFiles_='/store/user/johnda/AODSIM/Py8PtGun_bb_noPU_AODSIM/180731_210318/0000/step_AODSIM_1.root'
#inputFiles_='/store/user/johnda/AODSIM/Py8PtGun_bb_noPU_AODSIM/180731_210318/0000/step_AODSIM_2.root'
#inputFiles_='file:/eos/uscms/store/group/lpcml/QCDToGG_Pt_80_120_13TeV_TuneCUETP8M1_noPU_AODSIM/180809_215549/0000/step_full_1.root'
#inputFiles_='root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRFall17DR/RSGluonToTT_M-4000_TuneCUETP8M1_14TeV-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/0021E007-AEB9-E711-9AD1-E0071B7AC700.root'
#inputFiles_='root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/30002/FEB0DA54-FDB2-E711-AE70-0242AC110002.root'
inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/eusai/step2_ttbar_p8_03/step2_qcd8_1109.root'
#inputFiles_='file:../step_full.root'

#maxEvents_=100
#skipEvents_=0#
outputFile_ = 'test.root'

cmd="cmsRun %s inputFiles=%s outputFile=%s" %(cfg,inputFiles_,outputFile_)
print '%s'%cmd
os.system(cmd)
