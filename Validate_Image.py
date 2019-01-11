import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

f = h5py.File('output/test/QCD.hdf5')
print(list(f.keys()))

idx = 9
img_crop = f['X_jets'][idx,...]

plt.imshow(img_crop[...,2], vmin=1.e-3, cmap='Greys', norm=LogNorm(), alpha=0.9) # Tracks
plt.imshow(img_crop[...,1], vmin=1.e-3, cmap='Blues', norm=LogNorm(), alpha=0.9) # ECAL
plt.imshow(img_crop[...,0], vmin=1.e-3, cmap='Oranges', norm=LogNorm(), alpha=0.9) #HCAL

# axis labeling
# Note: due to the way imshow() renders images by default, the below lines will appear to 'flip' the image
ax = plt.axes()
plt.xlim([0., 125.+1.])
plt.xlabel(r'$\mathrm{i\varphi}$', size=14)
ax.xaxis.set_tick_params(direction='in', which='major', length=6.)
plt.ylim([0., 125.+1.])
plt.ylabel(r'$\mathrm{i\eta}$', size=14)
ax.yaxis.set_tick_params(direction='in', which='major', length=6.)
plt.show()

'''
ecal_showers = f['X_jets'][:,...,1]
plt.imshow(np.mean(ecal_showers, axis=0), vmin=1.e-3, vmax=ecal_showers.max(), cmap='hot_r', norm=LogNorm())
plt.colorbar(fraction=0.046, pad=0.04)
plt.show()
'''
