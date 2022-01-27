#!/usr/bin/env python3

##########################################################
# A collection of scripts for testing the behavior of ISC.
#
# Each script (segments between the large comment headers)
# can be run independently after the imports.
#
# In each script, ISC is tested in terms of how accurately it
# recovers the original coefficients.
#
# Each script deals with slightly different scenarios
# (testing on continuous data, testing on concatenated
# segments taken from a larger set, etc).
#
##########################################################


import numpy as np
import sys
import matplotlib.pyplot as plt
# module from hypercodes
sys.path.append('/home/adamb/hypercodes/current_code/')
import couplingISC as isc



#########################################################
# Baseline ISC behavior as a function of SNR
#
# ISC is tested here on a long continuous data set,
# with varying SNR levels
#########################################################

# params for random data generation
tr_no = 750
vox_no = 1000
data_dims = (tr_no, vox_no)
shift = 3
noise_range = [0.1, 10, 0.25]
noise_levels = np.arange(noise_range[0], noise_range[1], noise_range[2])
beta = np.asarray([0, 0.2, 0, 0.2, 0.5, 2, 0])  # make sure the number of coeffs in beta is shift*2+1
beta = beta/np.sum(np.absolute(beta))  # normalize coefficients
# rep_no = 10

# preallocate accuracy var
acc = np.zeros((len(noise_levels), 2))

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
f = plt.errorbar(np.log10(snr), acc[:, 0], yerr=acc[:, 1], fmt='o-')
plt.title('Coeff detection error as a function of SNR')
plt.xlabel('log10 SNR')
plt.ylabel('Median summed error of coeffs')
plt.savefig('baseline_ISC.png')



############################################################
# Version with non-continuous, concatenated segments
#
# ISC is tested here on concatenated segments from a long,
# continuous data set, as a function of SNR.
############################################################

# params for random data generation
tr_no = 1500
vox_no = 1000
data_dims = (tr_no, vox_no)
shift = 3
noise_range = [0.1, 10, 0.25]
noise_levels = np.arange(noise_range[0], noise_range[1], noise_range[2])
beta = np.asarray([0, 0.2, 0, 0.2, 0.5, 2, 0])  # make sure the number of coeffs in beta is shift*2+1
beta = beta/np.sum(np.absolute(beta))  # normalize coefficients

# params for segmentation
bounds = np.arange(0, 1500, 50)
bounds1 = bounds[0::2]
bounds2 = bounds[1::2]

# preallocate accuracy var
acc = np.zeros((len(noise_levels), 2))

# iterate over noise levels
noise_it = 0
for noise_level in noise_levels:

    print('Noise level: ' + str(noise_level))

    # generate random data
    X, Y = isc.get_coupled_set_2d(data_dims, beta, shift, noise_level)

    # get segments
    for i in range(len(bounds1)):
        if i == 0:
            X1 = X[bounds1[i]:bounds2[i], :]
            Y1 = Y[bounds1[i]:bounds2[i], :]
        else:
            X1 = np.concatenate((X1, X[bounds1[i]:bounds2[i], :]), axis=0)
            Y1 = np.concatenate((Y1, X[bounds1[i]:bounds2[i], :]), axis=0)

    # calculate ISC
    beta_est = isc.coupling_onestep(X1, Y1, shift)

    # get difference in betas (~ accuracy)
    acc[noise_it, 0] = np.median(np.sqrt(np.sum((beta_est - beta)**2, 1)))  # median srt of sum of squared diffs
    acc[noise_it, 1] = np.std(np.sqrt(np.sum((beta_est - beta) ** 2, 1)))  # std

    # increment
    noise_it = noise_it + 1


# plot the results
snr = 1/noise_levels
f = plt.errorbar(np.log10(snr), acc[:, 0], yerr=acc[:, 1], fmt='o-')
plt.title('Coeff detection error as a function of SNR, concatenated segments!')
plt.xlabel('log10 SNR')
plt.ylabel('Median summed error of coeffs')
plt.savefig('concatenated_segments_ISC.png')



############################################################
# Version with non-continuous segments, as with speech turns.
#
# ISC is tested here on segments from a long, continuous
# data set. The segments are NOT concatenated, instead, ISC
# is estimated on each segment and the average over their
# results is used for final accuracy.
#
# ISC accuracy is reported / tested as a function of SNR.
############################################################

# params for random data generation
tr_no = 1500
vox_no = 1000
data_dims = (tr_no, vox_no)
shift = 3
noise_range = [0.1, 10, 0.25]
noise_levels = np.arange(noise_range[0], noise_range[1], noise_range[2])
beta = np.asarray([-0.3, 0.2, 0, 0.2, 0.5, 2, 0.7])  # make sure the number of coeffs in beta is shift*2+1
beta = beta/np.sum(np.absolute(beta))  # normalize coefficients

# params for segmentation
bounds = np.arange(0, 1500, 50)
bounds1 = bounds[0::2]
bounds2 = bounds[1::2]

# preallocate accuracy var
acc = np.zeros((len(noise_levels), len(bounds1), 2))

# iterate over noise levels
noise_it = 0
for noise_level in noise_levels:

    print('Noise level: ' + str(noise_level))

    # generate random data
    X, Y = isc.get_coupled_set_2d(data_dims, beta, shift, noise_level)

    # iterate over segments
    for i in range(len(bounds1)):
        X1 = X[bounds1[i]:bounds2[i], :]
        Y1 = Y[bounds1[i]:bounds2[i], :]

        # calculate ISC
        beta_est = isc.coupling_onestep(X1, Y1, shift)

        # get difference in betas (~ accuracy)
        acc[noise_it, i, 0] = np.median(np.sqrt(np.sum((beta_est - beta)**2, 1)))  # median srt of sum of squared diffs
        acc[noise_it, i, 1] = np.std(np.sqrt(np.sum((beta_est - beta) ** 2, 1)))  # std

    # increment
    noise_it = noise_it + 1

# average over segments
acc = np.mean(acc, 1)

# plot the results
snr = 1/noise_levels
f = plt.errorbar(np.log10(snr), acc[:, 0], yerr=acc[:, 1], fmt='o-')
plt.title('Coeff detection error as a function of SNR, averaged ISC from segments!')
plt.xlabel('log10 SNR')
plt.ylabel('Median summed error of coeffs')
plt.savefig('averaged_over_segments_ISC.png')




############################################################
# Version with non-continuous segments, as with speech turns.
# BUT with fixed SNR, testing across different segment lengths (TRs)
#
# ISC is tested here on segments from a long, continuous
# data set. The segments are NOT concatenated, instead, ISC
# is estimated on each segment and the average over their
# results is used for final accuracy.
#
# ISC accuracy is reported / tested as a function of TRs
# in each segment.
############################################################

# params for random data generation
tr_no = 1500
vox_no = 1000
data_dims = (tr_no, vox_no)
shift = 3
noise_level = 1
beta = np.asarray([-0.3, 0.2, 0, 0.2, 0.5, 2, 0.7])  # make sure the number of coeffs in beta is shift*2+1
beta = beta/np.sum(np.absolute(beta))  # normalize coefficients
segment_lengths = np.arange(10, 50, 2)

# params for segmentation
bounds = np.arange(0, 1500, np.max(segment_lengths))
bounds1 = bounds[0::2]

# preallocate accuracy var
acc = np.zeros((len(segment_lengths), len(bounds1), 2))

# generate random data
X, Y = isc.get_coupled_set_2d(data_dims, beta, shift, noise_level)

# iterate over segment length values
segment_tr_it = 0
for segment_length in segment_lengths:

    print('Segment length (TRs): ' + str(segment_length))

    # iterate over segments
    for i in range(len(bounds1)):
        X1 = X[bounds1[i]:bounds1[i]+segment_length, :]
        Y1 = Y[bounds1[i]:bounds1[i]+segment_length, :]

        # calculate ISC
        beta_est = isc.coupling_onestep(X1, Y1, shift)

        # get difference in betas (~ accuracy)
        acc[segment_tr_it, i, 0] = np.median(np.sqrt(np.sum((beta_est - beta)**2, 1)))  # median srt of sum of squared diffs
        acc[segment_tr_it, i, 1] = np.std(np.sqrt(np.sum((beta_est - beta) ** 2, 1)))  # std

    # increment
    segment_tr_it = segment_tr_it + 1

# average over segments
acc = np.mean(acc, 1)

# plot the results
f = plt.errorbar(segment_lengths, acc[:, 0], yerr=acc[:, 1], fmt='o-')
plt.title('Coeff detection error as a function of TRs, averaged ISC from segments, SNR=' + str(1/noise_level))
plt.xlabel('Segment length (TRs)')
plt.ylabel('Median summed error of coeffs')
plt.savefig('ISC_per_trs.png')










