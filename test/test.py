import numpy as np
import healpy as hp
import sys
sys.path.insert(0, '../build')
import healpeak as hpk

nside = 128
npix = hp.nside2npix(nside)

m = np.random.standard_normal([1000, npix])

bins = np.linspace(-3, 3, 20)

print(hpk.peak_counts(m, bins))
print(hpk.void_counts(m, bins))