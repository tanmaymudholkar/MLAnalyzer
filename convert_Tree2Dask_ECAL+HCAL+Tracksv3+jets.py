import numpy as np
import ROOT
from root_numpy import tree2array
from dask.delayed import delayed
import dask.array as da
from skimage.measure import block_reduce
import os

#eosDir='/eos/uscms/store/user/mba2012/IMGs'
#eosDir='~/work/MLHEP/CMSSW_8_0_26_patch1/src/ggAnalysis/ggNtuplizer/test'
#eosDir='/eos/cms/store/user/mandrews/OPENDATA/IMGs/MGG90_Eta23'
eosDir='/uscms/home/bburkle/nobackup/working_area/CMSSW_9_3_0/src/MLAnalyzer/output/test'
outDir='/uscms/home/bburkle/nobackup/working_area/CMSSW_9_3_0/src/MLAnalyzer/output/test'
#decays = ['DiPhotonBorn_MGG90_Eta23', 'GluGluHToGG_MGG90_Eta23', 'GJet_MGG90_Eta23']
decays = ['TTbar','QCD']
#decays = ['dummy']

#chunk_size_ = 250
chunk_size_ = 100
#scale = [100., 150.]
scale = [1., 1.]
jet_shape = 125

@delayed
def load_X(tree, start_, stop_, branches_, readouts, scale):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    # Convert the object array X to a multidim array:
    # 1: for each event x in X, concatenate the object columns (branches) into a flat array of shape (readouts*branches)
    # 2: reshape the flat array into a stacked array: (branches, readouts)
    # 3: embed each stacked array as a single row entry in a list via list comprehension
    # 4: convert this list into an array with shape (events, branches, readouts) 
    X = np.array([np.concatenate(x).reshape(len(branches_),readouts[0]*readouts[1]) for x in X])
    #print "X.shape:",X.shape
    X = X.reshape((-1,len(branches_),readouts[0],readouts[1]))
    X = np.transpose(X, [0,2,3,1])

    # Rescale
    X /= scale 
    return X

@delayed
def load_single(tree, start_, stop_, branches_):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_)
    X = np.array([x[0] for x in X])
    return X

@delayed
def load_vector(tree, start_, stop_, branches_, idx_=None):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_)
    X = np.array([np.concatenate(x) for x in X])
    if idx_ is not None:
        X = X[:,idx_]
    return X

@delayed
def load_X_upsampled(tree, start_, stop_, branches_, readouts, scale, upscale):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    # Convert the object array X to a multidim array:
    # 1: for each event x in X, concatenate the object columns (branches) into a flat array of shape (readouts*branches)
    # 2: reshape the flat array into a stacked array: (branches, readouts)
    # 3: embed each stacked array as a single row entry in a list via list comprehension
    # 4: convert this list into an array with shape (events, branches, readouts) 
    X = np.array([np.concatenate(x).reshape(len(branches_),readouts[0]*readouts[1]) for x in X])
    #print "X.shape:",X.shape
    X = X.reshape((-1,len(branches_),readouts[0],readouts[1]))

    #print "unsampled.shape",X.shape
    X = np.stack([tile_stacked_array(x, upscale) for x in X])
    #print "upsampled.shape",X.shape
    X = np.transpose(X, [0,2,3,1])

    # Rescale
    X /= scale 
    return X

from numpy.lib.stride_tricks import as_strided

def tile_stacked_array(X, upscale):
    #print "un-tile_stacked.shape",X.shape
    X = np.stack([tile_array(x, upscale, upscale) for x in X])
    #print "tile_stacked.shape",X.shape
    return X
    
def tile_array(x, b0, b1):
    r, c = x.shape                                    # number of rows/columns
    rs, cs = x.strides                                # row/column strides 
    x = as_strided(x, (r, b0, c, b1), (rs, 0, cs, 0)) # view a as larger 4D array
    return x.reshape(r*b0, c*b1)                      # create new 2D array
#    return x.reshape(r*b0, c*b1)/(b0*b1)              # create new 2D array

def block_resample_EE(X):

    return np.array([resample_EE(x) for x in X])

def resample_EE(imgECAL, factor=2):
    
    imgECAL = np.squeeze(imgECAL)
    #print('imgECAL.shape:',imgECAL.shape)
    
    # EE-
    imgEEm = imgECAL[:140-85] # EE- in the first 55 rows
    imgEEm = np.pad(imgEEm, ((1,0),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEm_dn = block_reduce(imgEEm, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEm_dn_up = tile_array(imgEEm_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor**2
    imgECAL[:140-85] = imgEEm_dn_up[1:] ## replace the old EE- rows
    
    # EE+
    imgEEp = imgECAL[140+85:] # EE+ in the last 55 rows
    imgEEp = np.pad(imgEEp, ((0,1),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEp_dn = block_reduce(imgEEp, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEp_dn_up = tile_array(imgEEp_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor*factor
    imgECAL[140+85:] = imgEEp_dn_up[:-1] # replace the old EE+ rows
    
    return np.expand_dims(imgECAL, -1)

def crop_jet(imgECAL, iphi, ieta):
    iphi = int(iphi*5 + 2)
    ieta = int(ieta*5 + 2)
    off = jet_shape//2
    #print('iphi:%d, ieta:%d'%(iphi, ieta))
    if iphi < off:
        #print('1')
        diff = off-iphi
        #diff = diff+150
        #print(diff)
        img_crop = np.concatenate((imgECAL[ieta-off:ieta+off+1,-diff:],
                                   imgECAL[ieta-off:ieta+off+1,:iphi+off+1]), axis=1)
    elif 360-iphi < off:
        #print('2')
        diff = off - (360-iphi)
        #diff = diff+150
        #print('diff:',diff)
        img_crop = np.concatenate((imgECAL[ieta-off:ieta+off+1,iphi-off:], 
                                   imgECAL[ieta-off:ieta+off+1,:diff+1]), axis=1)
    else:
        #print('0')
        img_crop = imgECAL[ieta-off:ieta+off+1,iphi-off:iphi+off+1]
    #print(img_crop.shape)
    return img_crop

@delayed
def crop_jet_block(Xs, iphis, ietas):
    X = np.array([crop_jet(x,iphi,ieta) for x,iphi,ieta in zip(Xs,iphis,ietas)])
    #print X.shape
    return X

for j,decay in enumerate(decays):

    if j == 0:
    #if j == 0 or j == 1:
        pass
        #continue

    #tfile_str = 'output_dijet.root'
    tfile_str = '%s/%s_test.root'%(eosDir,decay)
    #tfile_str = 'output_ggqq.root'
    #tfile_str = 'output_numEvent20.root'
    #tfile_str = '%s/%s_IMG.root'%(eosDir,decay)
    #tfile_str = '%s/ggtree_mc_single.root'%(eosDir)
    #tfile_str = '%s/ggtree_mc.root'%(eosDir)
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    #tree = tfile.Get('ggNtuplizer/EventTree')
    nevts = tree.GetEntries()
    #neff = (nevts//1000)*1000
    #neff = (nevts//100)*100
    #neff = 84600
    #neff = 175000
    #neff = 135600
    #neff = 200
    neff = int(nevts)
    chunk_size = chunk_size_
    if neff < chunk_size:
      chunk_size = neff
    if neff > nevts:
        neff = int(nevts)
        chunk_size = int(nevts)
    if neff > chunk_size:
        assert neff%chunk_size == 0
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    ## eventId
    ##branches = ["event"]
    #branches = ["eventId"]
    #eventId = da.concatenate([\
    #            da.from_delayed(\
    #                load_single(tree,i,i+chunk_size, branches),\
    #                shape=(chunk_size,),\
    #                dtype=np.int32)\
    #            for i in range(0,neff,chunk_size)])
    #print " >> %s: %s"%(branches[0],eventId.shape)

    ## runId
    #branches = ["runId"]
    #runId = da.concatenate([\
    #            da.from_delayed(\
    #                load_single(tree,i,i+chunk_size, branches),\
    #                shape=(chunk_size,),\
    #                dtype=np.int32)\
    #            for i in range(0,neff,chunk_size)])
    #print " >> %s: %s"%(branches[0],runId.shape)

    # ECAL
    readouts = [280,360]
    branches = ["ECAL_energy"]
    X_ECAL = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[0]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s: %s"%(branches[0],X_ECAL.shape)

    # ECAL with resampled EE
    X_ECAL_EEup = X_ECAL.map_blocks(lambda x: block_resample_EE(x), dtype=np.float32)
    print " >> %s: %s"%('ECAL_EEup_energy',X_ECAL_EEup.shape)

    # Tracks at ECAL
    readouts = [280,360]
    branches = ["ECAL_tracksPt"]
    #branches = ["ECAL_tracksQPt"] # for Qxpt weighted 
    #branches = ["ECAL_EndtracksPt"] # for pt weighted at ECAL face
    X_TracksAtECAL = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[0]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s: %s"%(branches[0],X_ECAL.shape)

    # HBHE upsample
    readouts = [56,72]
    branches = ["HBHE_energy"]
    upscale = 5
    X_HBHE_up = da.concatenate([\
                da.from_delayed(\
                    load_X_upsampled(tree,i,i+chunk_size, branches, readouts, scale[1], upscale),\
                    shape=(chunk_size, readouts[0]*upscale, readouts[1]*upscale, len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s(upsampled): %s"%(branches[0],X_HBHE_up.shape)
    
    #X_MuonsAtECAL, 
    X_ECAL_stacked = da.concatenate([X_TracksAtECAL, X_ECAL_EEup, X_HBHE_up], axis=-1)
    print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked.shape)

    # EB
    readouts = [170,360]
    #branches = ["HBHE_energy_EB"]
    #branches = ["TracksQPt_EB","EB_energy"]
    #branches = ["TracksPt_EB","EB_energy"]
    branches = ["EndTracksPt_EB","EB_energy"]
    #branches = ["EB_energy"]
    #branches = ["EB_energy","HBHE_energy_EB","Tracks_EB"]
    X_EB = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[0]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s: %s"%(branches[0],X_EB.shape)

    # EE-
    readouts = [100,100]
    branches = ["TracksPt_EEm","EEm_energy","HBHE_energy_EEm"]
    #branches = ["EndTracksPt_EEm","EEm_energy","HBHE_energy_EEm"] # for pt weighted at ECAL face
    #branches = ["TracksQPt_EEm","EEm_energy","HBHE_energy_EEm"] # for Qxpt weighted 
    #branches = ["EEm_energy","HBHE_energy_EEm","Tracks_EEm"]
    X_EEm = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[1]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s: %s"%(branches[0],X_EEm.shape)

    # EE+
    readouts = [100,100]
    branches = ["TracksPt_EEp","EEp_energy","HBHE_energy_EEp"]
    #branches = ["EndTracksPt_EEp","EEp_energy","HBHE_energy_EEp"] # for pt weighted at ECAL face
    #branches = ["TracksQPt_EEp","EEp_energy","HBHE_energy_EEp"] # for Qxpt weighted
    #branches = ["EEp_energy","HBHE_energy_EEp","Tracks_EEp"]
    X_EEp = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[1]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s: %s"%(branches[0],X_EEp.shape)

    # HBHE
    readouts = [56,72]
    branches = ["HBHE_energy"]
    X_HBHE = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[1]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s: %s"%(branches[0],X_HBHE.shape)

    # HBHE_EM
    #readouts = [56,72]
    #branches = ["HBHE_EMenergy"]
    #X_HBHE_EM = da.concatenate([\
    #            da.from_delayed(\
    #                load_X(tree,i,i+chunk_size, branches, readouts, scale[1]),\
    #                shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
    #                dtype=np.float32)\
    #            for i in range(0,neff,chunk_size)])
    #print " >> %s: %s"%(branches[0],X_HBHE_EM.shape)

    # HB_EB upsample
    readouts = [34,72]
    branches = ["HBHE_energy_EB"]
    upscale = 5
    X_HBHE_EB_up = da.concatenate([\
                da.from_delayed(\
                    load_X_upsampled(tree,i,i+chunk_size, branches, readouts, scale[1], upscale),\
                    shape=(chunk_size, readouts[0]*upscale, readouts[1]*upscale, len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s(upsampled): %s"%(branches[0],X_HBHE_EB_up.shape)

    X_EB = da.concatenate([X_EB, X_HBHE_EB_up], axis=-1)
    print " >> %s: %s"%('X_EB', X_EB.shape)

    for ijet in [0]: #[0,1]:

        print ' >> jet index',ijet

        # jet m0
        branches = ["jetM"]
        jetM = da.concatenate([\
                    da.from_delayed(\
                        load_vector(tree,i,i+chunk_size,branches,ijet),\
                        shape=(chunk_size,),\
                        dtype=np.float32)\
                    for i in range(0,neff,chunk_size)])
        print "  >> jetM:", jetM.shape

        # jet pt
        branches = ["jetPt"]
        jetPt = da.concatenate([\
                    da.from_delayed(\
                        load_vector(tree,i,i+chunk_size,branches,ijet),\
                        shape=(chunk_size,),\
                        dtype=np.float32)\
                    for i in range(0,neff,chunk_size)])
        print "  >> jetPt:", jetPt.shape

        # jet seed iphi
        branches = ["jetSeed_iphi"]
        jetSeed_iphi = da.concatenate([\
                    da.from_delayed(\
                        load_vector(tree,i,i+chunk_size,branches,ijet),\
                        shape=(chunk_size,),\
                        dtype=np.float32)\
                    for i in range(0,neff,chunk_size)])
        print "  >> jetSeed_iphi:", jetSeed_iphi.shape

        # jet seed ieta
        branches = ["jetSeed_ieta"]
        jetSeed_ieta = da.concatenate([\
                    da.from_delayed(\
                        load_vector(tree,i,i+chunk_size,branches,ijet),\
                        shape=(chunk_size,),\
                        dtype=np.float32)\
                    for i in range(0,neff,chunk_size)])
        print "  >> jetSeed_ieta:", jetSeed_ieta.shape

        # jet window
        X_jets = da.concatenate([\
                    da.from_delayed(\
                        crop_jet_block(X_ECAL_stacked[i:i+chunk_size], jetSeed_iphi[i:i+chunk_size], jetSeed_ieta[i:i+chunk_size]),\
                        shape=(chunk_size, jet_shape, jet_shape, 3),\
                        dtype=np.float32)\
                    for i in range(0,neff,chunk_size)])

        # Class label
        label = j
        #label = 1
        print "  >> Class label:",label
        y = da.from_array(\
                #np.full(len(eventId), label, dtype=np.float32),\
                np.full(len(X_jets), label, dtype=np.float32),\
                chunks=(chunk_size,))
        print "  >> y shape:",y.shape

        #file_out_str = "%s/test_dijet_ijet%d.hdf5"%(outDir,ijet)
        file_out_str = "%s/%s.hdf5"%(outDir,decay)
        #file_out_str = "test%d_numEvent1.hdf5"%label
        #file_out_str = "%s/%s_IMGall_RH%d_n%d_label%d.hdf5"%(eosDir,decay,int(scale[0]),neff,label)
        #file_out_str = "%s/%s_IMG_RH%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),neff//1000.)
        #file_out_str = "%s/%s_IMG_EBEEHBup_RH%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),neff//1000.)
        #file_out_str = "%s/%s_IMG_RH%d-%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),int(scale[1]),neff//1000.)
        print "  >> Writing to:", file_out_str
        if os.path.isfile(file_out_str):
            os.remove(file_out_str)
        #da.to_hdf5(file_out_str, {'/X_EB': X_EB, 'X_EEm': X_EEm, 'X_EEp': X_EEp, 'X_HBHE': X_HBHE, '/y': y}, compression='lzf')
        #da.to_hdf5(file_out_str, {'/X': X_EB, 'X_EEm': X_EEm, 'X_EEp': X_EEp, 'X_HBHE': X_HBHE, '/y': y}, compression='lzf')
        da.to_hdf5(file_out_str, {
                                  #'eventId': eventId,
                                  #'runId': runId,
                                  'X_ECAL': X_ECAL,
                                  #'X_ECAL_EEup': X_ECAL_EEup,
                                  'X_ECAL_stacked': X_ECAL_stacked,
                                  'X_EB': X_EB,
                                  'X_EEm': X_EEm,
                                  'X_EEp': X_EEp,
                                  'X_HBHE': X_HBHE,
                                  #'X_HBHE_EM': X_HBHE_EM,
                                  'X_HBHE_EB_up': X_HBHE_EB_up,
                                  'jetSeed_iphi': jetSeed_iphi,
                                  'jetSeed_ieta': jetSeed_ieta,
                                  'X_jets': X_jets,
                                  'jetPt': jetPt,
                                  'jetM': jetM,
                                  '/y': y
                                  }, compression='lzf')

        print "  >> Done.\n"
