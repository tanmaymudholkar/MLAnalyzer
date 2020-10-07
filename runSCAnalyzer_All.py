import os
from glob import glob
import re
import argparse

parser = argparse.ArgumentParser(description='Run SCAnalyzer')
parser.add_argument('-d','--decay', required=True, help='Decay:Single*Pt50',type=str)
args = parser.parse_args()

eosDir='/eos/uscms/store/user/mba2012'
decay='%s_AODSIM'%args.decay

cfg='RecHitAnalyzer/python/SCAnalyzer_cfg.py'
inputFiles_ = ['file:%s'%path for path in glob('%s/AODSIM/%s/*/*/step*root'%(eosDir,decay))]

listname = 'list_%s.txt'%decay
with open(listname, 'w') as list_file:
    for inputFile in inputFiles_:
        list_file.write("%s\n" % inputFile)

maxEvents_=-1
#maxEvents_=5000
skipEvents_=0

cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=%s/IMGs/%s_IMG.root"%(cfg,listname,maxEvents_,skipEvents_,eosDir,decay)
#cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=test_IMG.root"%(cfg,listname,maxEvents_,skipEvents_)
#print '%s'%cmd
os.system(cmd)
#test
#os.system('mv cEB*.eps %s/'%(inputTag))
