import h5py
import numpy

f = h5py.File("test/conversion_test.hdf5", 'r')
print f.keys()
n=f.get('jetPt')
print n
jetPt=numpy.array(n)
print jetPt
print jetPt.shape
print numpy.array(f['jet_ngen'])
print numpy.array(f['jet_npf'])
print numpy.array(f['jet_genpart_pdgid'])
print numpy.array(f['jet_genpart_charge'])