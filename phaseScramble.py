import numpy as np
import matplotlib.pyplot as mplot
from numpy.random import default_rng
from math import pi


def get_random_matrix(d0, d1):
    '''
    Get a 2D array filled with random values between -1 and 1
    Params "d0" and "d1" specify the size of the array.
    '''
    rng = default_rng()
    rand_data = rng.random((d0, d1))*2-1
    return rand_data

# get a random array in range -1, 1
voxelNo = 50000
trNo = 1000
data = get_random_matrix(trNo, voxelNo)*2-1

# use fft.rfft as our data is real
data_rfft = np.fft.rfft(data, axis=0)

# do fft.irfft
data_irfft = np.fft.irfft(data_rfft, axis=0)

# check for numerical inaccuracies
diffs = (np.abs(data-data_irfft)).flatten()
print('Maximum difference after forward-backward DFT pass: ' + '{:.3e}'.format(diffs.max()))

# check covariance
datacov = np.cov(data, rowvar=False)

# random phases
rng = default_rng()
rand_phase = (rng.random((trNo, 1))*2*pi-pi)*1j

# # histogram of differences
# mplot.hist(diffs, bins=100)



