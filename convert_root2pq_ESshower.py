import pyarrow.parquet as pq
import pyarrow as pa # pip install pyarrow==0.7.1
import ROOT
import numpy as np
from numpy.lib.stride_tricks import as_strided
import glob, os

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('-i', '--infile', default='output.root', type=str, help='Input root file.')
parser.add_argument('-i', '--infile', default=['output.root'], nargs='+', type=str, help='Input root file.')
#parser.add_argument('-i', '--infile', default=['output.root'], type=list, help='Input root file.')
#parser.add_argument('-o', '--outdir', default='parquet', type=str, help='Output pq file dir.')
parser.add_argument('-o', '--outdir', default='.', type=str, help='Output pq file dir.')
parser.add_argument('-d', '--decay', default='test', type=str, help='Decay name.')
parser.add_argument('-n', '--idx', default=0, type=int, help='Input root file index.')
#parser.add_argument('-w', '--wgt_file', default=None, type=str, help='Weight file.')
parser.add_argument('-w', '--wgt_files', default=None, nargs='+', type=str, help='Weight file.')
args = parser.parse_args()

nXY = 40
nSTRIP = 32

def tile_array(x, b0, b1):
    r, c = x.shape                                    # number of rows/columns
    rs, cs = x.strides                                # row/column strides
    x = as_strided(x, (r, b0, c, b1), (rs, 0, cs, 0)) # view a as larger 4D array
    return x.reshape(r*b0, c*b1)/(b0*b1)              # create new 2D array

def tile_ES(es_img, plane):
    assert (plane == 'X') or (plane == 'Y')
    if plane == 'X':
        es_img = np.float32(es_img).reshape(nXY, nXY*nSTRIP)
        es_img = tile_array(es_img, nSTRIP, 1)
    else:
        es_img = np.float32(es_img).reshape(nXY*nSTRIP, nXY)
        es_img = tile_array(es_img, 1, nSTRIP)
    return es_img

def crop_ESshower(es_img, ix, iy, window=7*nSTRIP):
    off = window//2
    seed_x = int(ix)*nSTRIP + nSTRIP//2 # upper side of center for even window
    seed_y = int(iy)*nSTRIP + nSTRIP//2 # upper side of center for even window
    #es_shower = es_img[seed_row-window-nSTRIP:seed_row+window,seed_col-window-nSTRIP:seed_col+window]
    es_shower = es_img[seed_y-off:seed_y+off,seed_x-off:seed_x+off] # x:cols, y:rows
    return es_shower

def crop_EEshower(ee_img, ix, iy, window=16):
    off = window//2
    seed_x = int(ix) # upper side of center for even window
    seed_y = int(iy) # upper side of center for even window
    ee_shower = ee_img[seed_y-off:seed_y+off,seed_x-off:seed_x+off] # x:cols, y:rows
    return ee_shower

def pa_array(d):
    arr = pa.array([d]) if np.isscalar(d) or type(d) == list else pa.array([d.tolist()])
    #print(arr.type)
    ## double to single float
    if arr.type == pa.float64():
        arr = arr.cast(pa.float32())
    elif arr.type == pa.list_(pa.float64()):
        arr = arr.cast(pa.list_(pa.float32()))
    elif arr.type == pa.list_(pa.list_(pa.float64())):
        arr = arr.cast(pa.list_(pa.list_(pa.float32())))
    elif arr.type == pa.list_(pa.list_(pa.list_(pa.float64()))):
        arr = arr.cast(pa.list_(pa.list_(pa.list_(pa.float32()))))
    elif arr.type == pa.list_(pa.list_(pa.list_(pa.list_(pa.float64())))):
        arr = arr.cast(pa.list_(pa.list_(pa.list_(pa.list_(pa.float32())))))
    #else:
    #    print('Unknown type for conversion to (list of) floats',arr.type)
    #print(arr.type)
    return arr

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

rhTreeStr = args.infile
print " >> Input file:",rhTreeStr
rhTree = ROOT.TChain("fevt/RHTree")
for f in rhTreeStr:
  rhTree.Add(f)
nEvts = rhTree.GetEntries()
assert nEvts > 0
print " >> nEvts:",nEvts
outStr = '%s/%s.parquet.%d'%(args.outdir, args.decay, args.idx)
print " >> Output file:",outStr

##### EVENT SELECTION START #####

# Event range to process
iEvtStart = 0
#iEvtEnd   = 50
iEvtEnd   = nEvts
assert iEvtEnd <= nEvts
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

nPhos = 0
data = {} # Arrays to be written to parquet should be saved to data dict
es_planes = {}
es_shower = {}
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

    # Only select events with same number of photon candidates identified in EE and ES
    # 1.479 < |eta_EE| < 3.0
    # 1.653 < |eta_ES| < 2.6
    if len(rhTree.seed_ix) != len(rhTree.SC_iz): continue
    # Check EE
    nPhoEE = 0
    for i in range(len(rhTree.SC_iz)):
        ix, iy, iz = rhTree.SC_iphi[i], rhTree.SC_ieta[i], rhTree.SC_iz[i]
        # Boundary cuts to ensure photon seed remains on-center
        if ( ix < 8 or ix+7 > 99 ): continue
        if ( iy < 8 or iy+7 > 99 ): continue
        nPhoEE += 1
    # Check ES
    nPhoES = 0
    for i in range(len(rhTree.seed_ix)):
        ix, iy, iz = rhTree.seed_ix[i], rhTree.seed_iy[i], rhTree.seed_iz[i]
        # Boundary cuts to ensure photon seed remains on-center
        if ( ix < 3 or ix > 40-4 ): continue
        if ( iy < 3 or iy > 40-4 ): continue
        nPhoES += 1
    if nPhoEE != nPhoES: continue

    nPhoEvt = nPhoEE
    if args.wgt_files is not None:
        rands = np.random.random((nPhoEvt, nPasses))

    for i in range(nPhoEvt):

        data['idx'] = idx + [i]
        #data['m0'] = rhTree.m0

        data['m'] = rhTree.SC_mass[i]
        data['pt'] = rhTree.SC_pT[i]
        data['iphi'] = rhTree.SC_iphi[i]
        data['ieta'] = rhTree.SC_ieta[i]
        data['crystal_maxE_X'] = rhTree.SC_X[i]
        data['crystal_maxE_Y'] = rhTree.SC_Y[i]
        data['iz'] = rhTree.SC_iz[i]
        data['SC_genZ'] = rhTree.SC_genZ[i]
        data['SC_daughter1_projEE'] = [rhTree.SC_daughter1_projEE_X[i], rhTree.SC_daughter1_projEE_Y[i]]
        data['SC_daughter1_pT'] = rhTree.SC_daughter1_pT[i]
        data['SC_daughter2_projEE'] = [rhTree.SC_daughter2_projEE_X[i], rhTree.SC_daughter2_projEE_Y[i]]
        data['SC_daughter2_pT'] = rhTree.SC_daughter2_pT[i]

        #if rhTree.pho_pT[i] > 100.: continue

        if args.wgt_files is not None:
            keepEG = True
            for p in range(nPasses):
                if rands[i,p] < get_weight_2d(data['m'], data['pt'], m_edgess[p], pt_edgess[p], hmvpts[p]):
                    keepEG = False
            if keepEG == False: continue

        data['A_p4'] = [
            rhTree.SC_E[i]
            ,rhTree.SC_pT[i]
            ,rhTree.SC_eta[i]
            ,rhTree.SC_phi[i]
            ]

        data['A_projEE'] = [rhTree.SC_projEE_X[i], rhTree.SC_projEE_Y[i]]
        data['A_pT'] = rhTree.SC_pT[i]

        data['pho_p4'] = [rhTree.pho_E[i], rhTree.pho_pT[i], rhTree.pho_eta[i], rhTree.pho_phi[i]]
        data['pho_id'] = [
                rhTree.pho_r9[i]
                ,rhTree.pho_sieie[i]
                ,rhTree.pho_phoIso[i]
                ,rhTree.pho_chgIso[i]
                ,rhTree.pho_chgIsoWrongVtx[i]
                ,rhTree.pho_Eraw[i]
                ,rhTree.pho_phiWidth[i]
                ,rhTree.pho_etaWidth[i]
                ,rhTree.pho_scEta[i]
                ,rhTree.pho_sieip[i]
                ,rhTree.pho_s4[i]
            ]
        #data['pho_vars'] = [
        #        rhTree.pho_r9[i]
        #        ,rhTree.pho_HoE[i]
        #        ,rhTree.pho_hasPxlSeed[i]
        #        ,rhTree.pho_sieie[i]
        #        ,rhTree.pho_phoIso[i]
        #        ,rhTree.pho_trkIso[i]
        #        ,rhTree.pho_chgIsoCorr[i]
        #        ,rhTree.pho_neuIsoCorr[i]
        #        ,rhTree.pho_phoIsoCorr[i]
        #        ,rhTree.pho_bdt[i]
        #    ]

        ##### Make shower images

        # EE energy at EE geometry
        ix, iy, iz = rhTree.SC_iphi[i], rhTree.SC_ieta[i], rhTree.SC_iz[i]
        ee_img = rhTree.EEp_energy if iz > 0 else rhTree.EEm_energy
        ee_img = np.float32(ee_img).reshape(100,100)
        data['X_ee'] = crop_EEshower(ee_img, ix, iy).reshape(1,16,16)

        # ES geometry images

        # ES energy at ES geometry
        ix, iy, iz = rhTree.seed_ix[i], rhTree.seed_iy[i], rhTree.seed_iz[i]
        # If photon candidate is on positive endcap
        if iz > 0:
            es_planes['X'] = rhTree.ESpX_energy
            es_planes['Y'] = rhTree.ESpY_energy
        # If photon candidate is on negative endcap
        else:
            es_planes['X'] = rhTree.ESmX_energy
            es_planes['Y'] = rhTree.ESmY_energy

        # Make image for each ES plane of the appropriate endcap
        for k in ['X', 'Y']:
            # NOTE: Arrays are stored by strip!
            # Since each ES plane has strips segmented only along a single direction,
            # input arrays will have total length = n_sensors*n_strips_per_sensor * n_sensors.
            # These must be reshaped differently for the X:F, Y:R planes
            # then only upsampled along the unsegmented direction
            # to output the same 1280 x 1280 square image
            es_img = tile_ES(es_planes['%s'%(k)], k)
            es_shower[k] = crop_ESshower(es_img, ix, iy).reshape(1,224,224)

        # EE energy at ES geometry
        eeAtes_img = rhTree.EE_energy_ESp if iz > 0 else rhTree.EE_energy_ESm
        eeAtes_img = np.float32(eeAtes_img).reshape(1280,1280)
        eeAtes_shower = crop_ESshower(eeAtes_img, ix, iy).reshape(1,224,224)

        # Tracks pt at ES geometry
        #tkAtes_img = rhTree.TracksPt_ESp if iz > 0 else rhTree.TracksPt_ESm
        #tkAtes_img = np.float32(tkAtes_img).reshape(1280,1280)
        #tkAtes_shower = crop_ESshower(tkAtes_img, ix, iy).reshape(1,224,224)

        data['X_cms'] = np.concatenate([es_shower['X'], es_shower['Y'], eeAtes_shower], axis=0)
        #print(data['X_cms'].shape)

        #pqdata = [pa.array([d]) if np.isscalar(d) or type(d) == list else pa.array([d.tolist()]) for d in data.values()]
        pqdata = [pa_array(d) for d in data.values()]
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
X = pqIn.read_row_group(0, columns=['idx.list.item','iphi','ieta','iz','pho_p4.list.item']).to_pydict()
X = pqIn.read_row_group(1, columns=['idx.list.item','iphi','ieta','iz','pho_p4.list.item']).to_pydict()
print(X)
#X = pqIn.read_row_group(0, columns=['X.list.item.list.item.list.item']).to_pydict()['X']
#X = pqIn.read(['X.list.item.list.item.list.item']).to_pydict()['X']
#X = np.float32(X)
