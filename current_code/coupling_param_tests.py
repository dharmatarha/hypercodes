#!/usr/bin/env python3

import numpy as np
import sys
import matplotlib.pyplot as plt
# module from hypercodes
sys.path.append('/home/adamb/hypercodes/current_code/')
import couplingISC as isc

# params for random data generation
tr_no = 1500
vox_no = 1000
data_dims = (tr_no, vox_no)
shift = 3
# noise_level = 0.1
noise_range = [0, 5, 0.1]
noise_levels = np.arange(noise_range[0], noise_range[1], noise_range[2])
beta = np.asarray([0, 0.2, 0, 0.2, 0.5, 2, 0])
beta = beta/np.sum(beta)  # normalize coefficients
# rep_no = 10

# preallocate accuracy var
acc = np.zeros((50, 2))

# iterate over noise levels
noise_it = 0
for noise_level in noise_levels:

    print('Noise level: ' + str(noise_level))

#    # iterate over repetitions
#    for rep in range(rep_no):

    # generate random data
    X, Y = isc.get_coupled_set_2d(data_dims, beta, shift, noise_level)

    # calculate ISC
    beta_est = isc.coupling_onestep(X, Y, shift)

    # get difference in betas (~ accuracy)
    acc[noise_it, 0] = np.median(np.sqrt(np.sum((beta_est - beta)**2, 1)))  # median srt of sum of squared diffs
    acc[noise_it, 1] = np.std(np.sqrt(np.sum((beta_est - beta) ** 2, 1))) # std

    # increment
    noise_it = noise_it + 1

# plot the results
snr = 1/noise_levels
f = plt.errorbar(snr, acc[:, 0], yerr=acc[:, 1], fmt='o-')
plt.title('Coeff detection error as a function of SNR')
plt.xlabel('SNR')
plt.ylabel('Median summed error of coeffs')