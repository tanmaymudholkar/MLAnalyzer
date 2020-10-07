import os
from glob import glob
import re
import argparse

parser = argparse.ArgumentParser(description='Run SCRegressor')
#parser.add_argument('-d','--decay', required=True, help='Decay:Single*Pt50',type=str)
parser.add_argument('-m','--mass', default='Bkg',type=str)
parser.add_argument('-s','--sel', required=True,type=str)
args = parser.parse_args()

#eosDir='/eos/uscms/store/user/mba2012'
#eosDir='/eos/uscms/store/user/lpcml/mandrews'
#xrootd='root://cmsxrootd.fnal.gov' # FNAL
xrootd='root://cms-xrd-global.cern.ch' #CERN
#xrootd='root://eoscms.cern.ch' # CERN
#decay='%s'%args.decay
decay='HAHMHToAA_AToGG'
#decay='SUSYGluGluToHToAA_AToGG'
mass='%s'%args.mass
sel='%s'%args.sel
#mass='10GeV'
#evtcont = 'MINIAODSIM'
evtcont = 'AODSIM'

cfg='RecHitAnalyzer/python/SCRegressor_cfg.py'
#inputFiles_ = ['%s/%s'%(xrootd,path) for path in glob( '/store/mc/RunIIFall17MiniAODv2/SUSYGluGluToHToAA_AToGG_M-15_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/*root')]#645076A5-B87D-E811-8CEE-A0369FD0B342.root'
#inputFiles_ = ['file:%s'%path for path in glob('%s/%s/%s/*/*/step*root'%(eosDir, evtcont, decay))]
#inputFiles_ = ['%s/%s'%(xrootd,path) for path in glob('%s/AODSIM/%s/*/*/step*root'%(eosDir,decay))]
#inputFiles_ = ['%s/%s'%(xrootd,path.replace('/eos/uscms','')) for path in glob('%s/%s/%s/*/*/step*root'%(eosDir, evtcont, decay))]
#print(inputFiles_)

listname = 'LISTS/list_Ma-%s_%s_%s.txt'%(mass,evtcont, decay) #signal
#listname='LISTS/list_DiphotonJets_MGG-80toInf_%s.txt' %sel #bkg
#with open(listname, 'w') as list_file:
#    for inputFile in inputFiles_:
#        list_file.write("%s\n" % inputFile)

maxEvents_=2000
#maxEvents_=100000
skipEvents_=0

#decay=decay.replace('_%s'%evtcont,'')
#subdirout = 'CUTS_GEN'
#subdirout = 'CUTS_KIN'
#subdirout = 'CUTS_LOOSE'
#subdirout = 'dPhidEta'
#subdirout = 'Pi0_flat_mvpt'
cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=Output/%s_M%s_output.root" %(cfg,listname,maxEvents_,skipEvents_,sel,mass)
#cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=test_IMG.root"%(cfg,listname,maxEvents_,skipEvents_)
#print '%s'%cmd
os.system(cmd)

#os.system('mv cEB*.eps %s/'%(inputTag))
