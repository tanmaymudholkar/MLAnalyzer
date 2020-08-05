import ROOT
import numpy as np

tree = ROOT.TChain("fevt/RHTree")
#tree.Add("output_OctoPi0_e60_m200_ctau2em5To12_eta0To1p4_noPU_AODSIM.root")
tree.Add("output.root")
#tree.Draw('SC_dvtx', 'nPhoEvent >= 2')
#tree.Draw('SC_dvtx')
#tree.Draw('nRecHits:SC_dvtx')
#tree.Draw('nPhoEvent')
#tree.Draw('SC_dvtx:nPhoEvent')
#tree.Draw('SC_recoIdx')
#tree.Draw('SC_dvtx:pho_E')
#tree.Draw('SC_dvtx:nPhoEvent', 'nPhoEvent >= 2')
#tree.Draw('pho_eta:nPhoEvent', 'nPhoEvent >= 2')
#tree.Draw('pho_eta:SC_dvtx', 'nPhoEvent >= 2')
#tree.Draw('pho_r9:SC_dvtx', 'nPhoEvent >= 2')
brs = list(tree.GetListOfBranches())
brs = [br.GetName() for br in brs]
print(brs)

h, c = {}, {}
wd, ht = 680, 480

k = 'nPhoEvent'
h[k] = ROOT.TH1F(k, k, 9, 0, 9.)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'SC_dvtxGen'
h[k] = ROOT.TH1F(k, k, 70, 0, 350)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'SC_dvtx'
h[k] = ROOT.TH1F(k, k, 50, 0, 250)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'SC_vtxT'
h[k] = ROOT.TH1F(k, k, 50, 0, 250)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'SC_vtxZ'
h[k] = ROOT.TH1F(k, k, 50, 0, 250)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'SC_boost'
h[k] = ROOT.TH1F(k, k, 50, 250., 350.)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'pho_r9'
h[k] = ROOT.TH1F(k, k, 48, 0, 1.2)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'pho_sieie'
h[k] = ROOT.TH1F(k, k, 50, 0, 0.05)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'eb_time'
h[k] = ROOT.TH1F(k, k, 50, -25., 25.)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'SC_vtxTvvtxZ'
h[k] = ROOT.TH2F(k, k, 10, 0, 250, 14, 0., 350.)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'SC_vtxTvvtxZGen'
h[k] = ROOT.TH2F(k, k, 10, 0, 250, 14, 0., 350.)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'r9vdvtx'
h[k] = ROOT.TH2F(k, k, 10, 0, 250, 24, 0., 1.2)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'ptvdvtx'
h[k] = ROOT.TH2F(k, k, 10, 0, 250, 50, 0., 100.)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'etavdvtx'
h[k] = ROOT.TH2F(k, k, 10, 0, 250, 28, -1.4, 1.4)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'dRvdvtx'
h[k] = ROOT.TH2F(k, k, 10, 0, 250, 20, 0., 5.)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'timevdvtx'
h[k] = ROOT.TH2F(k, k, 10, 0, 250, 25, -25., 25.)
c[k] = ROOT.TCanvas(k, k, wd, ht)
k = 'nPhoEventvdvtxGen'
h[k] = ROOT.TH2F(k, k, 10, 0, 250, 9, 0., 9.)
c[k] = ROOT.TCanvas(k, k, wd, ht)

print(tree.GetEntries())
nEvts = tree.GetEntries()
for i in range(nEvts):
    tree.GetEntry(i)

    nAGenValid = len(tree.SC_dvtx)
    if nAGenValid == 0: continue

    h['nPhoEvent'].Fill(tree.nPhoEvent)
    #print('tree.SC_dvtx:%d, tree.SC_phi:%d, tree.SC_iphi:%d, nPhoEvent:%d, nPhoRegress:%d'\
    #    %(len(tree.SC_dvtx), len(tree.SC_phi), len(tree.SC_iphi), tree.nPhoEvent, tree.nPhoRegress))
    #print(list(tree.SC_dvtx))
    h['nPhoEventvdvtxGen'].Fill(tree.SC_dvtx[0], tree.nPhoEvent)
    #eb = np.float32(tree.EB_time)
    #eb = eb[eb != 0.]
    #for t in eb:
    #    #if t == 0.: continue
    #    h['eb_time'].Fill(t)
    #break
    for a in range(nAGenValid):
        h['SC_vtxTvvtxZGen'].Fill(tree.SC_vtxZ[a], tree.SC_vtxT[a])
        h['SC_dvtxGen'].Fill(tree.SC_dvtx[a])
        if tree.SC_recoIdx[a] == -1: continue
        h['SC_vtxTvvtxZ'].Fill(tree.SC_vtxZ[a], tree.SC_vtxT[a])
        h['SC_dvtx'].Fill(tree.SC_dvtx[a])
        h['SC_vtxT'].Fill(tree.SC_vtxT[a])
        h['SC_vtxZ'].Fill(tree.SC_vtxZ[a])
        h['SC_boost'].Fill(tree.SC_boost[a])
        reco_idx = int(tree.SC_recoIdx[a])
        h['r9vdvtx'].Fill(tree.SC_dvtx[a], tree.pho_r9[reco_idx])
        h['ptvdvtx'].Fill(tree.SC_dvtx[a], tree.SC_pT[a])
        h['etavdvtx'].Fill(tree.SC_dvtx[a], tree.SC_eta[a])
        h['dRvdvtx'].Fill(tree.SC_dvtx[a], abs(tree.SC_DR[a])/0.0174)
        sc_time = np.float32(tree.SC_time)
        sc_time = sc_time[sc_time != 0.]
        for t in sc_time:
            h['eb_time'].Fill(t)
            h['timevdvtx'].Fill(tree.SC_dvtx[a], t)
    for p in range(tree.nPhoRegress):
        h['pho_r9'].Fill(tree.pho_r9[p])
        h['pho_sieie'].Fill(tree.pho_sieie[p])
    #if i >= 200: break

print(type(ROOT.TH2F()))

def normalizeTH2(k):
    nX = {}
    for ix in range(1, h[k].GetNbinsX()+1):
        nX[ix] = sum([h[k].GetBinContent(ix, iy_) for iy_ in range(1, h[k].GetNbinsY()+1)])

    for ix in range(1, h[k].GetNbinsX()+1):
        for iy in range(1, h[k].GetNbinsY()+1):
            binc = h[k].GetBinContent(ix, iy)
            binc = binc/nX[ix] if nX[ix] > 0. else binc
            h[k].SetBinContent(ix, iy, binc)

#h['nPhoEvent'].Draw()
for k in h.keys():
    print(k)
    c[k].cd()
    h[k].SetTitle(k)
    if 'time' in k:
        c[k].SetLogz()
    if type(h[k]) == type(ROOT.TH2F()) and 'vtxTvvtxZ' not in k:
        #normalizeTH2(k)
        h[k].Draw('COL Z')
    else:
        h[k].Draw()
    c[k].Draw()
    c[k].Update()
    #c[k].Print('Plots/normed/%s.eps'%k)
    #c[k].Print('Plots/nonnormed/%s.eps'%k)
