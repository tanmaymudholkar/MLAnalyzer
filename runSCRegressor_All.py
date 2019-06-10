import os
from glob import glob
import re
import argparse

parser = argparse.ArgumentParser(description='Run SCRegressor')
parser.add_argument('-d','--decay', required=True, help='Decay:Single*Pt50',type=str)
args = parser.parse_args()

#eosDir='/eos/uscms/store/user/mba2012'
eosDir='/eos/uscms/store/user/lpcml/mandrews'
xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
decay='%s'%args.decay

evtcont = 'MINIAODSIM'
#evtcont = 'AODSIM'

cfg='RecHitAnalyzer/python/SCRegressor_cfg.py'
inputFiles_ = ['file:%s'%path for path in glob('%s/%s/%s/*/*/step*root'%(eosDir, evtcont, decay))]
#inputFiles_ = ['%s/%s'%(xrootd,path) for path in glob('%s/AODSIM/%s/*/*/step*root'%(eosDir,decay))]
#inputFiles_ = ['%s/%s'%(xrootd,path.replace('/eos/uscms','')) for path in glob('%s/%s/%s/*/*/step*root'%(eosDir, evtcont, decay))]
#print(inputFiles_)

listname = 'LISTS/list_%s_%s.txt'%(evtcont, decay)
with open(listname, 'w') as list_file:
    for inputFile in inputFiles_:
        list_file.write("%s\n" % inputFile)

maxEvents_=-1
maxEvents_=100000
skipEvents_=0

decay=decay.replace('_%s'%evtcont,'')
#cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=%s/IMG/CUTS_GEN/%s_IMG.root"%(cfg,listname,maxEvents_,skipEvents_,eosDir,decay)
#cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=%s/IMG/CUTS_KIN/%s_IMG.root"%(cfg,listname,maxEvents_,skipEvents_,eosDir,decay)
cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=%s/IMG/dPhidEta/%s_IMG.root"%(cfg,listname,maxEvents_,skipEvents_,eosDir,decay)
#cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=test_IMG.root"%(cfg,listname,maxEvents_,skipEvents_)
#print '%s'%cmd
os.system(cmd)

#os.system('mv cEB*.eps %s/'%(inputTag))
