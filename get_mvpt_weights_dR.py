import ROOT
import numpy as np
np.random.seed(0)
import glob, os

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-d', '--genDR', default=10, type=int, help='gen-level dR.')
parser.add_argument('-p', '--p_drop', default=1., type=float, help='p(drop) scale.')
args = parser.parse_args()

genDR = args.genDR
p_drop_scale = args.p_drop

xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
eosDir='/eos/uscms/store/user/lpcml/mandrews/IMG'

#pu = 'noPU'
#pu = '2016_25ns_Moriond17MC_PoissonOOTPU'
#decay = 'DoublePi0Pt15To100_m0To1600_pythia8_%s_mlog_ptexp'%pu
#decay = '%s_genDR%d_recoDR16_seedPos_phoVars_IMG'%(decay, genDR)
#date_str = '190310_173649' # DR10
pu = 'PU2017'
decay = 'DoublePi0Pt10To100_m0To1600_pythia8_PU2017_MINIAODSIM_ext2'
date_str = '190522_041152'

# Paths to input files
rhFileList = '%s/%s/%s/*/output_*.root'%(eosDir, decay, date_str)
print(" >> Input file list: %s"%rhFileList)
rhFileList = glob.glob(rhFileList)
assert len(rhFileList) > 0
print(" >> %d files found"%len(rhFileList))
rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList]
print(' >> Input File[0]: %s'%rhFileList[0])

rhTree = ROOT.TChain("fevt/RHTree")
for f in rhFileList:
    rhTree.Add(f)
nEvts = rhTree.GetEntries()
#rhTree.SetBranchStatus("EB_*", 0)
#rhTree.SetBranchStatus("SC_energy*", 0)
#rhTree.SetBranchStatus("SC_time", 0)
assert nEvts > 0
print " >> nEvts:",nEvts

def get_weight_2d(m0, pt, m0_edges, pt_edges, wgts):
    idx_m0 = np.argmax(m0 <= m0_edges)-1
    idx_pt = np.argmax(pt <= pt_edges)-1
    #print(idx_m0, idx_pt)
    return wgts[idx_m0, idx_pt]

##### EVENT SELECTION START #####

nPasses = 2

# Event range to process
iEvtStart = 0
iEvtEnd   = 50000
iEvtEnd   = nEvts
assert iEvtEnd <= nEvts
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

eb_xtal = 0.0174

sw = ROOT.TStopwatch()
sw.Start()
d = {}
hmvpts, m_edgess, pt_edgess = {}, {}, {}
for p in range(nPasses):

    print " >> Pass:", p
    nAcc = 0
    sc_mass_, sc_pT_, sc_dR_ = [], [], []
    hnPho = ROOT.TH2F("hnPho%d"%p, "hnPho;m_{#pi^{0}};p_{T,#pi^{0}}", 16, 0, 1.6 , 20, 20., 100.)

    for iEvt in range(iEvtStart,iEvtEnd):

        rhTree.GetEntry(iEvt)

        if iEvt % 10000 == 0:
            print " .. Processing entry",iEvt

        #d['pho_pt']          = list(rhTree.pho_pT)
        #d['pho_r9']          = list(rhTree.pho_r9)
        #d['pho_sieie']       = list(rhTree.pho_sieie)
        #d['pho_phoIso']      = list(rhTree.pho_phoIso)
        #d['pho_ecalIsoCorr'] = list(rhTree.pho_ecalIsoCorr)
        #d['pho_trkIso']      = list(rhTree.pho_trkIso)
        #d['pho_HoE']         = list(rhTree.pho_HoE)
        #d['pho_hasPxlSeed']  = list(rhTree.pho_hasPxlSeed)

        d['SC_mass'] = list(rhTree.SC_mass)
        d['SC_pT'] = list(rhTree.SC_pT)
        d['SC_dR'] = list(rhTree.SC_DR)

        nPhoEvt = len(d['SC_mass'])
        if p > 0:
            rands = np.random.random(nPhoEvt)

        for i in range(nPhoEvt):

            if d['SC_pT'][i] < 20.: continue
            if d['SC_dR'][i]/0.0174 > 10.: continue

            #if d['pho_r9'][i] <= 0.5: continue
            #if d['pho_HoE'][i] >= 0.08: continue
            #if d['pho_r9'][i] <= 0.85:
            #    if d['pho_sieie'][i] >= 0.015: continue
            #    if d['pho_phoIso'][i] >= 4.: continue
            #    if d['pho_trkIso'][i] >= 6.: continue
            #if d['pho_hasPxlSeed'][i] == True: continue

            sc_mass = d['SC_mass'][i]
            sc_pT = d['SC_pT'][i]
            sc_dR = d['SC_dR'][i]

            if p > 0:
                keepEG = True
                for p_ in range(p):
                    if rands[i] < get_weight_2d(sc_mass, sc_pT, m_edgess[p_], pt_edgess[p_], hmvpts[p_]):
                        keepEG = False
                if keepEG == False: continue

            hnPho.Fill(sc_mass, sc_pT)
            nAcc += 1

            sc_mass_.append(sc_mass)
            sc_pT_.append(sc_pT)
            sc_dR_.append(sc_dR)

    sw.Stop()
    print " >> Real time:",sw.RealTime()/60.,"minutes"
    print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
    print " >> nPi0s: ", nAcc

    sc_mass_ = np.array(sc_mass_)
    sc_pT_ = np.array(sc_pT_)
    sc_dR_ = np.array(sc_dR_)

    hmvpt, m_edges, pt_edges = np.histogram2d(sc_mass_, sc_pT_, range=((0., 1.6), (20.,100.)), bins=(16, 20))
    hmvpt = hmvpt/hmvpt.max()
    floor = np.mean(hmvpt.flatten()) - 2.*np.std(hmvpt.flatten())
    if floor < 0.:
        print " >> Forcing floor to %f -> 0."%floor
        floor = 0.
    print " >> w(m_v_pt) min, floor:",hmvpt.min(),floor
    hmvpt[hmvpt < floor] = floor
    hmvpt = hmvpt/hmvpt.max()
    hmvpt = hmvpt-hmvpt.min()
    print " >> p(drop) scale:",p_drop_scale
    hmvpt = p_drop_scale*hmvpt
    print " >> w(m_v_pt)*p(drop) min, max:",hmvpt.min(), hmvpt.max()
    wgt_file = 'WEIGHTS/%s_mvpt_weights_pdrop%.2f_pass%d.npz'%(decay, p_drop_scale, p)
    np.savez(wgt_file,\
            mvpt = hmvpt,
            m_edges = m_edges,
            pt_edges = pt_edges,
            sc_mass = sc_mass_,
            sc_pT = sc_pT_,
            sc_dR = sc_dR_
            )
    w = np.load(wgt_file)
    hmvpts[p], m_edgess[p], pt_edgess[p] = w['mvpt'], w['m_edges'], w['pt_edges']

    #cnPhow = ROOT.TCanvas("cnPhow%d"%p, "cnPhow%d"%p, 600, 600)
    #hnPho.Draw("COL Z")

hnPhow = ROOT.TH2F("hnPhow", "hnPhow;m_{#pi^{0}};p_{T,#pi^{0}}", 16, 0, 1.6 , 20, 20., 100.)
nAcc = 0
for iEvt in range(iEvtStart,iEvtEnd):

    rhTree.GetEntry(iEvt)

    if iEvt % 10000 == 0:
        print " .. Processing entry",iEvt

    #d['pho_pt']          = list(rhTree.pho_pT)
    #d['pho_r9']          = list(rhTree.pho_r9)
    #d['pho_sieie']       = list(rhTree.pho_sieie)
    #d['pho_phoIso']      = list(rhTree.pho_phoIso)
    #d['pho_ecalIsoCorr'] = list(rhTree.pho_ecalIsoCorr)
    #d['pho_trkIso']      = list(rhTree.pho_trkIso)
    #d['pho_HoE']         = list(rhTree.pho_HoE)
    #d['pho_hasPxlSeed']  = list(rhTree.pho_hasPxlSeed)

    d['SC_mass'] = list(rhTree.SC_mass)
    d['SC_pT'] = list(rhTree.SC_pT)
    d['SC_dR'] = list(rhTree.SC_DR)

    nPhoEvt = len(d['SC_mass'])
    #rands = np.random.random(nPhoEvt)
    rands = np.random.random((nPhoEvt, nPasses))

    for i in range(nPhoEvt):

        if d['SC_pT'][i] < 20.: continue
        if d['SC_dR'][i]/0.0174 > 10.: continue

        #if d['pho_r9'][i] <= 0.5: continue
        #if d['pho_HoE'][i] >= 0.08: continue
        #if d['pho_r9'][i] <= 0.85:
        #    if d['pho_sieie'][i] >= 0.015: continue
        #    if d['pho_phoIso'][i] >= 4.: continue
        #    if d['pho_trkIso'][i] >= 6.: continue
        #if d['pho_hasPxlSeed'][i] == True: continue

        sc_mass = d['SC_mass'][i]
        sc_pT = d['SC_pT'][i]
        #sc_dR = d['SC_dR'][i]

        keepEG = True
        for p in range(nPasses):
            if rands[i,p] < get_weight_2d(sc_mass, sc_pT, m_edgess[p], pt_edgess[p], hmvpts[p]):
                keepEG = False
        #if rands[i] < get_weight_2d(sc_mass, sc_pT, m_edges, pt_edges, hmvpt):
        #    continue
        if keepEG == False: continue

        hnPhow.Fill(sc_mass, sc_pT)
        nAcc += 1

print " >> nPi0s: ", nAcc
#cnPhow = ROOT.TCanvas("cnPhow", "cnPhow", 600, 600)
#hnPhow.Draw("COL Z")

hFile = ROOT.TFile("WEIGHTS/%s_hnPho_pdrop%.2f_passes%d.root"%(decay, p_drop_scale, nPasses),"RECREATE")

hnPho.SetMinimum(0.)
hnPho.Write()
hnPhoX = hnPho.ProjectionX()
hnPhoX.GetYaxis().SetRangeUser(0., 1.2*hnPhoX.GetMaximum())
hnPhoX.Write()
hnPhoY = hnPho.ProjectionY()
hnPhoY.GetYaxis().SetRangeUser(0., 1.2*hnPhoY.GetMaximum())
hnPhoY.Write()

hnPhow.SetMinimum(0.)
hnPhow.Write()
hnPhowX = hnPhow.ProjectionX()
hnPhowX.GetYaxis().SetRangeUser(0., 1.2*hnPhowX.GetMaximum())
hnPhowX.Write()
hnPhowY = hnPhow.ProjectionY()
hnPhowY.GetYaxis().SetRangeUser(0., 1.2*hnPhowY.GetMaximum())
hnPhowY.Write()

hFile.Close()
