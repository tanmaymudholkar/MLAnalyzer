import pyarrow.parquet as pq
import pyarrow as pa # pip install pyarrow==0.7.1
import ROOT
import numpy as np
np.random.seed(0)
import glob, os

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('-i', '--infile', default='output.root', type=str, help='Input root file.')
parser.add_argument('-i', '--infile', default=['output.root'], nargs='+', type=str, help='Input root file.')
#parser.add_argument('-i', '--infile', default=['output.root'], type=list, help='Input root file.')
parser.add_argument('-o', '--outdir', default='.', type=str, help='Output pq file dir.')
parser.add_argument('-d', '--decay', default='test', type=str, help='Decay name.')
parser.add_argument('-n', '--idx', default=0, type=int, help='Input root file index.')
#parser.add_argument('-w', '--wgt_file', default=None, type=str, help='Weight file.')
parser.add_argument('-w', '--wgt_files', default=None, nargs='+', type=str, help='Weight file.')
args = parser.parse_args()

def crop_EBshower(imgEB, iphi, ieta, window=32):

    # NOTE: image window here should correspond to the one used in RHAnalyzer
    off = window//2
    iphi = int(iphi)+1 # seed positioned at [15,15]
    ieta = int(ieta)+1 # seed positioned at [15,15]

    # Wrap-around on left side
    if iphi < off:
        diff = off-iphi
        img_crop = np.concatenate((imgEB[:,ieta-off:ieta+off,-diff:],
                                   imgEB[:,ieta-off:ieta+off,:iphi+off]), axis=-1)
    # Wrap-around on right side
    elif 360-iphi < off:
        diff = off - (360-iphi)
        img_crop = np.concatenate((imgEB[:,ieta-off:ieta+off,iphi-off:],
                                   imgEB[:,ieta-off:ieta+off,:diff]), axis=-1)
    # Nominal case
    else:
        img_crop = imgEB[:,ieta-off:ieta+off,iphi-off:iphi+off]

    return img_crop

def get_weight_2d(m0, pt, m0_edges, pt_edges, wgts):
    idx_m0 = np.argmax(m0 <= m0_edges)-1
    idx_pt = np.argmax(pt <= pt_edges)-1
    #print(idx_m0, idx_pt)
    return wgts[idx_m0, idx_pt]

hmvpts, m_edgess, pt_edgess = {}, {}, {}
if args.wgt_files is not None:
    nPasses = len(args.wgt_files)
    for p,wgt_file in enumerate(args.wgt_files):
        w = np.load(wgt_file)
        hmvpts[p], m_edgess[p], pt_edgess[p] = w['mvpt'], w['m_edges'], w['pt_edges']
    #w = np.load(args.wgt_file)
    #hmvpt, m_edges, pt_edges = w['mvpt'], w['m_edges'], w['pt_edges']

rhTreeStr = args.infile
print " >> Input file:",rhTreeStr
rhTree = ROOT.TChain("fevt/RHTree")
for f in rhTreeStr:
  rhTree.Add(f)
nEvts = rhTree.GetEntries()
assert nEvts > 0
print " >> nEvts:",nEvts
outStr = '%s/%s.parquet.%d'%(args.outdir, args.decay, args.idx)
#outStr = '%s/%s.reg_2reco.parquet.%d'%(args.outdir, args.decay, args.idx)
print " >> Output file:",outStr

##### EVENT SELECTION START #####

# Event range to process
iEvtStart = 0
iEvtEnd   = 10000
iEvtEnd   = nEvts
assert iEvtEnd <= nEvts
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

nPhos = 0
data = {} # Arrays to be written to parquet should be saved to data dict
sw = ROOT.TStopwatch()
sw.Start()
for iEvt in range(iEvtStart,iEvtEnd):

    # Initialize event
    rhTree.GetEntry(iEvt)

    if iEvt % 10000 == 0:
        print " .. Processing entry",iEvt

    #if rhTree.m0 < 100. or rhTree.m0 > 110.:
    #  continue

    idx = [rhTree.runId, rhTree.lumiId, rhTree.eventId]

    SC_energyT = rhTree.SC_energyT
    SC_energyZ = rhTree.SC_energyZ
    SC_energy  = rhTree.SC_energy

    pi0_mass = rhTree.SC_mass
    pi0_iphi = rhTree.SC_iphi
    pi0_ieta = rhTree.SC_ieta
    pi0_dR = rhTree.SC_DR

    pi0_E   = rhTree.SC_E
    pi0_pt  = rhTree.SC_pT
    pi0_phi = rhTree.SC_phi
    pi0_eta = rhTree.SC_eta

    pho_r9 =             rhTree.pho_r9
    pho_sieie =          rhTree.pho_sieie
    pho_phoIso =         rhTree.pho_phoIso
    pho_chgIso =         rhTree.pho_chgIso
    pho_chgIsoWrongVtx = rhTree.pho_chgIsoWrongVtx
    pho_Eraw =           rhTree.pho_Eraw
    pho_phiWidth =       rhTree.pho_phiWidth
    pho_etaWidth =       rhTree.pho_etaWidth
    pho_scEta =          rhTree.pho_scEta
    pho_sieip =          rhTree.pho_sieip
    pho_s4 =             rhTree.pho_s4

    pho_E   = rhTree.pho_E
    pho_pt  = rhTree.pho_pT
    pho_phi = rhTree.pho_phi
    pho_eta = rhTree.pho_eta

    #EB_time = np.array(rhTree.EB_time).reshape(170,360)
    #TracksAtEB_pt = np.array(rhTree.TracksPt_EB).reshape(170,360)
    #X_EB = np.stack([TracksAtEB_pt, EB_time], axis=0)
    #X_EB = np.array(rhTree.EB_energy).reshape(1,170,360)

    nPhoEvt = len(pi0_mass)
    rands = np.random.random((nPhoEvt, nPasses))
    for i in range(nPhoEvt):

        data['idx'] = idx + [i]
        #data['m0'] = rhTree.m0

        data['m'] = pi0_mass[i]
        data['pt'] = pi0_pt[i]
        data['iphi'] = pi0_iphi[i]
        data['ieta'] = pi0_ieta[i]
        data['dR'] = pi0_dR[i]

        #if data['ieta'] >= 170-16:
        #    continue

        if data['pt'] < 20.: continue
        if data['dR']/0.0174 > 10.: continue

        if args.wgt_files is not None:
            keepEG = True
            for p in range(nPasses):
                if rands[i,p] < get_weight_2d(data['m'], data['pt'], m_edgess[p], pt_edgess[p], hmvpts[p]):
                    keepEG = False
            if keepEG == False: continue
            #if rands[i] < get_weight_2d(data['m'], data['pt'], m_edges, pt_edges, hmvpt):
            #    continue

        data['pi0_p4'] = [pi0_E[i], pi0_pt[i], pi0_eta[i], pi0_phi[i]]
        data['pho_p4'] = [pho_E[i], pho_pt[i], pho_eta[i], pho_phi[i]]
        data['pho_id'] = [
            pho_r9[i]
            ,pho_sieie[i]
            ,pho_phoIso[i]
            ,pho_chgIso[i]
            ,pho_chgIsoWrongVtx[i]
            ,pho_Eraw[i]
            ,pho_phiWidth[i]
            ,pho_etaWidth[i]
            ,pho_scEta[i]
            ,pho_sieip[i]
            ,pho_s4[i]
        ]

        data['X'] = np.array(SC_energy[i]).reshape(1,32,32)
        sc_energyT = np.array(SC_energyT[i]).reshape(1,32,32)
        sc_energyZ = np.array(SC_energyT[i]).reshape(1,32,32)
        data['Xtz'] = np.concatenate((sc_energyT, sc_energyZ), axis=0)

        #sc_cms = crop_EBshower(X_EB, data['iphi'], data['ieta'])
        #if sc_cms.shape != data['Xtz'].shape:
        #    print(sc_cms.shape)
        #    print(data['Xtz'].shape)
        #data['Xcms'] = np.concatenate([sc_cms, data['Xtz']], axis=0)

        pqdata = [pa.array([d]) if np.isscalar(d) or type(d) == list else pa.array([d.tolist()]) for d in data.values()]
        table = pa.Table.from_arrays(pqdata, data.keys())

        if nPhos == 0:
            writer = pq.ParquetWriter(outStr, table.schema, compression='snappy')

        writer.write_table(table)

        nPhos += 1

writer.close()

sw.Stop()
print " >> nPhos:",nPhos
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
print " >> ======================================"

pqIn = pq.ParquetFile(outStr)
print(pqIn.metadata)
print(pqIn.schema)
#X = pqIn.read_row_group(0, columns=['m','pt','iphi','ieta','pt_reco']).to_pydict()
X = pqIn.read_row_group(0, columns=['idx.list.item','m','iphi','ieta','pi0_p4.list.item','pho_p4.list.item']).to_pydict()
print(X)
#X = pqIn.read_row_group(0, columns=['X.list.item.list.item.list.item']).to_pydict()['X']
#X = pqIn.read(['X.list.item.list.item.list.item']).to_pydict()['X']
#X = np.float32(X)
