import glob, os
import ROOT
import numpy as np
import h5py
from skimage.measure import block_reduce
from numpy.lib.stride_tricks import as_strided
import itertools

# TODO make it so we can read a file where each event has either 1 OR 2 jets, and process all jets

import argparse
parser = argparse.ArgumentParser(add_help=True, description='Process some integers.')
parser.add_argument('-l', '--label', default=1, type=int, help='Decay label.')
parser.add_argument('-n', '--fileNum', type=int, default=0, help='Choose which file you are processing')
parser.add_argument('-j', '--doJet', type=int, default=1, help='How many are in your root file. If 0, will do all jets. If 1 or 2, will only convert 1st or 2nd jet')
parser.add_argument('-g', '--granularity', default=1, type=int, help='Increased Image Pixel Granularity')
args = parser.parse_args()

def upsample_array(x, b0, b1):

    r, c = x.shape                                    # number of rows/columns
    rs, cs = x.strides                                # row/column strides
    x = as_strided(x, (r, b0, c, b1), (rs, 0, cs, 0)) # view as a larger 4D array

    return x.reshape(r*b0, c*b1)/(b0*b1)              # create new 2D array with same total occupancy 

def resample_EE(imgECAL, factor=2):

    # EE-
    imgEEm = imgECAL[:140-85] # EE- in the first 55 rows
    imgEEm = np.pad(imgEEm, ((1,0),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEm_dn = block_reduce(imgEEm, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEm_dn_up = upsample_array(imgEEm_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor**2
    imgECAL[:140-85] = imgEEm_dn_up[1:] ## replace the old EE- rows

    # EE+
    imgEEp = imgECAL[140+85:] # EE+ in the last 55 rows
    imgEEp = np.pad(imgEEp, ((0,1),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEp_dn = block_reduce(imgEEp, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEp_dn_up = upsample_array(imgEEp_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor*factor
    imgECAL[140+85:] = imgEEp_dn_up[:-1] # replace the old EE+ rows

    return imgECAL

def crop_jet(imgECAL, iphi, ieta, gran, jet_shape=125):

    # NOTE: jet_shape here should correspond to the one used in RHAnalyzer
    off = jet_shape//2
    iphi = int(iphi*5 + 2)*gran + gran//2 # 5 EB xtals per HB tower
    ieta = int(ieta*5 + 2)*gran + gran//2 # 5 EB xtals per HB tower

    # Wrap-around on left side
    if iphi < off:
        diff = off-iphi
        img_crop = np.concatenate((imgECAL[ieta-off:ieta+off+1,-diff:],
                                   imgECAL[ieta-off:ieta+off+1,:iphi+off+1]), axis=1)
    # Wrap-around on right side
    elif 360*gran-iphi < off:
        diff = off - (360*gran-iphi)
        img_crop = np.concatenate((imgECAL[ieta-off:ieta+off+1,iphi-off:],
                                   imgECAL[ieta-off:ieta+off+1,:diff+1]), axis=1)
    # Nominal case
    else:
        img_crop = imgECAL[ieta-off:ieta+off+1,iphi-off:iphi+off+1]

    return img_crop


outDir='/home/bjornb/fullSample_highGran/hdf5'
xrootd='root://cmsxrootd.fnal.gov' # FNAL

decays = ['QCD', 'TTbar']

scale = [1., 1.]
jet_shape = 125 * args.granularity
doJet = args.doJet
width = 280*args.granularity
height = 360*args.granularity

inputDir = '/home/bjornb/fullSample_highGran/ntuples'
if doJet == 1:
    ttbar_files = glob.glob('%s/ttbar/*' % inputDir)
    qcd_300to600 = glob.glob('%s/qcd_300to600/*' % inputDir)
    qcd_400to600 = glob.glob('%s/qcd_400to600/*' % inputDir)
    qcd_600to3000 = glob.glob('%s/qcd_600to3000/*' % inputDir)
elif doJet == 2:
    ttbar_files = glob.glob('%s/ttbar_2jets/*' % inputDir)
    qcd_300to600 = glob.glob('%s/qcd_300to600_2jets/*' % inputDir)
    qcd_400to600 = glob.glob('%s/qcd_400to600_2jets/*' % inputDir)
    qcd_600to3000 = glob.glob('%s/qcd_600to3000_2jets/*' % inputDir)

qcd_files = qcd_300to600 + qcd_400to600 + qcd_600to3000

TEST = True
test = 'test/conversion_test.hdf5'

# Loop over decays
for d, decay in enumerate(decays):

    if d != args.label:
        continue
   
    if d == 0:
        filelist = qcd_files
        tfile_idxs = range(1, len(qcd_files)+1)
    elif d == 1:
        filelist = ttbar_files
        tfile_idxs = range(1, len(ttbar_files)+1)
    else:
        print 'decay must be equal to 0 or 1'
        break
    	
    print '>> Doing decay[%d]: %s'%(d, decay)

    # Get root tree
    


    ############################################################################################################################
    ############################################################################################################################
    ################################O#V#E#R#R#I#D#E#############################################################################
    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################
    # tfile_str = str(filelist[args.fileNum])
    print " >> For input file:", tfile_str
    #tfile = ROOT.TXNetFile(tfile_str)
    tfile = ROOT.TFile('test/ttbar_new-production_test.root')
    # tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()

    if TEST:
        nevts = 100

    print " >> Total events:", nevts

    outPath = '%s/%s_IMGjet_%d-Granularity'%(outDir, decay, args.granularity)
    # note, if mutliple jobs are submitted to condor and multiple try to create the directory simulataneously, a lot of those jobs will fail
    if not os.path.isdir(outPath):
        os.makedirs(outPath)
    fout_str = '%s/%s_f%d_j%d_n%d.hdf5' % (outPath, decay, args.fileNum, doJet, nevts)
    if TEST:
        fout_str = test

    # create output file and define it's datasets
    fout = h5py.File(fout_str, 'w')
    fout.create_dataset('X_jets', (nevts, jet_shape, jet_shape, 8), compression='lzf') # note, number of image channels was hardcoded here
    fout.create_dataset('jetPt', (nevts, ), compression='lzf')
    fout.create_dataset('jetM', (nevts, ), compression='lzf')
    fout.create_dataset('y', (nevts, ), compression='lzf')

    # define branch names based on whether or not you are using high granularity images
    if args.granularity == 1:
        br_pt = 'ECAL_tracksPt_atECALfixIP'
        br_d0 = 'ECAL_tracksD0Sig_atECALfixIP'
        br_dz = 'ECAL_tracksDzSig_atECALfixIP'
        br_ecal = 'ECAL_energy'
        br_hcal = 'HBHE_energy'
        br_pix1 = 'BPIX_layer1_ECAL_atPV'
        br_pix2 = 'BPIX_layer2_ECAL_atPV'
        br_pix3 = 'BPIX_layer3_ECAL_atPV'
    else:
        br_pt = 'ECALadj_tracksPt_%dx%d'%(args.granularity, args.granularity)
        br_d0 = 'ECALadj_tracksD0Sig_%dx%d'%(args.granularity,args.granularity)
        br_dz = 'ECALadj_tracksDzSig_%dx%d'%(args.granularity,args.granularity)
        br_ecal = 'ECAL_energy'
        br_hcal = 'HBHE_energy'
        br_pix1 = 'BPIX_layer1_ECALadj_%dx%d'%(args.granularity,args.granularity)
        br_pix2 = 'BPIX_layer2_ECALadj_%dx%d'%(args.granularity,args.granularity)
        br_pix3 = 'BPIX_layer3_ECALadj_%dx%d'%(args.granularity,args.granularity)
    
    # convert to hdf5 file using an event loop
    for iEvt in range(0, nevts):
        
        if not iEvt % 100:
            print 'Processing event', iEvt
        if TEST:
            print iEvt
            
        # initialize tree
        tree.GetEntry(iEvt)

        # pyroot makes you normally do TreeName.BranchName
        # however, BranchName is a variable in our case, so we use
        # the getattr method
        TracksPt = np.array(getattr(tree, br_pt)).reshape(width, height)
        TracksD0 = np.array(getattr(tree, br_d0)).reshape(width, height)
        TracksDz = np.array(getattr(tree, br_dz)).reshape(width, height)
        Ecal = np.array(getattr(tree, br_ecal)).reshape(280,360)
        Ecal = resample_EE(Ecal)
        if args.granularity != 1:
            Ecal = upsample_array(Ecal, args.granularity, args.granularity)
        Hcal = np.array(getattr(tree, br_hcal)).reshape(56,72)
        Hcal = upsample_array(Hcal, 5*args.granularity, 5*args.granularity)
        pix1 = np.array(getattr(tree, br_pix1)).reshape(width, height)
        pix2 = np.array(getattr(tree, br_pix2)).reshape(width, height)
        pix3 = np.array(getattr(tree, br_pix3)).reshape(width, height)

        #jet_stack = np.stack((TracksPt, TracksD0, TracksDz, Ecal, Hcal, pix1, pix2, pix3), axis=-1)

        #del TracksPt
        #del TracksD0
        #del TracksDz
        #del Ecal
        #del Hcal
        #del pix1
        #del pix2
        #del pix3

        #TracksPt = None
        #TracksD0 = None
        #TracksDz = None
        #Ecal = None
        #Hcal = None
        #pix1 = None
        #pix2 = None
        #pix3 = None

        pts = tree.jetPt
        m0s = tree.jetM
        iphis = tree.jetSeed_iphi
        ietas = tree.jetSeed_ieta
        
        nJets = len(pts)

        for ijet in range(nJets):
            
            if doJet > 0 and ijet+1 != doJet:
                continue
            else:
                pass

            y = d
            pt = pts[ijet]
            m0 = m0s[ijet]
            iphi = iphis[ijet]
            ieta = ietas[ijet]
            
            # crop images individually so it is less cpu intensive
            TracksPt = crop_jet( TracksPt, iphi, ieta, args.granularity, jet_shape )
            TracksD0 = crop_jet( TracksD0, iphi, ieta, args.granularity, jet_shape )
            TracksDz = crop_jet( TracksDz, iphi, ieta, args.granularity, jet_shape )
            Ecal = crop_jet( Ecal, iphi, ieta, args.granularity, jet_shape )
            Hcal = crop_jet( Hcal, iphi, ieta, args.granularity, jet_shape )
            pix1 = crop_jet( pix1, iphi, ieta, args.granularity, jet_shape )
            pix2 = crop_jet( pix2, iphi, ieta, args.granularity, jet_shape )
            pix3 = crop_jet( pix3, iphi, ieta, args.granularity, jet_shape )
            
            X_jet = np.stack((TracksPt, TracksD0, TracksDz, Ecal, Hcal, pix1, pix2, pix3), axis=-1)

            TracksPt = None
            TracksD0 = None
            TracksDz = None
            Ecal = None
            Hcal = None
            pix1 = None
            pix2 = None
            pix3 = None

            # stacking images when we crop so we don't have two copies of all of the image channels saved as different variables
            #X_jet = crop_jet( np.concatenate([TracksPt, TracksD0, TracksDz, Ecal, Hcal, pix1, pix2, pix3], axis=-1), iphi, ieta, args.granularity )
            #X_jet = crop_jet( jet_stack, iphi, ieta, args.granularity )
            
            #del jet_stack
            #jet_stack = None

            fout['X_jets'][iEvt] = X_jet
            fout['jetPt'][iEvt] = pt
            fout['jetM'][iEvt] = m0
            fout['y'][iEvt] = y

            X_jets = None

        #TracksPt = None
        #TracksD0 = None
        #TracksDz = None
        #Ecal = None
        #Hcal = None
        #pix1 = None
        #pix2 = None
        #pix3 = None

    fout.close()

print "  >> Done.\n"