import os
from glob import glob
import re
import argparse

parser = argparse.ArgumentParser(description='Run SCRegressor')
parser.add_argument('-d','--decay', required=True, help='Decay:Single*Pt50',type=str)
args = parser.parse_args()

#eosDir='/eos/uscms/store/user/mba2012'
eosDir='/eos/uscms/store/user/lpcml/mandrews'
#xrootd='root://cmsxrootd.fnal.gov' # FNAL
xrootd='root://cmseos.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
decay='%s'%args.decay

#evtcont = 'MINIAODSIM'
evtcont = 'AODSIM'

#/eos/uscms/store/user/lpcml/mandrews/AODSIM/OctoPi0_e60_m200_ctau2em5To1p1e2_eta0To1p4_noPU_AODSIM/200728_145732/0000/step_AODSIM_noPU_ext4_1.root
cfg='RecHitAnalyzer/python/SCRegressor_cfg.py'
#inputFiles_ = ['file:%s'%path for path in glob('%s/%s/%s/*/*/step*root'%(eosDir, evtcont, decay))]
#inputFiles_ = ['%s/%s'%(xrootd,path) for path in glob('%s/AODSIM/%s/*/*/step*root'%(eosDir,decay))]
inputFiles_ = ['%s/%s'%(xrootd,path.replace('/eos/uscms','')) for path in glob('%s/%s/%s/*/*/step*root'%(eosDir, evtcont, decay))]
print(inputFiles_[0])

#'''
#listname = 'LISTS/list_%s_%s.txt'%(evtcont, decay)
listname = 'LISTS/list_%s.txt'%(decay)
with open(listname, 'w') as list_file:
    for inputFile in inputFiles_:
        list_file.write("%s\n" % inputFile)

maxEvents_=-1
maxEvents_=10
skipEvents_=0

#decay=decay.replace('_%s'%evtcont,'')
#subdirout = 'CUTS_GEN'
#subdirout = 'CUTS_KIN'
#subdirout = 'CUTS_LOOSE'
#subdirout = 'dPhidEta'
#subdirout = 'Pi0_flat_mvpt'
#cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=%s/IMG/%s/%s_IMG.root"%(cfg,listname,maxEvents_,skipEvents_,eosDir,subdirout,decay)
#cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=%s/IMG/%s/output.root"%(cfg,listname,maxEvents_,skipEvents_,eosDir,decay)
cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=output_%s.root"%(cfg,listname,maxEvents_,skipEvents_, decay)
#print '%s'%cmd
os.system(cmd)

#os.system('mv cEB*.eps %s/'%(inputTag))
#'''
