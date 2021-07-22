#!/usr/bin/env python3

'''
Script to test segmentation effects on ISC results.

The question was how to properly estimate ISC from segments of a continuous BOLD time series.
In our case, the continuous recording is from the Joint Storytelling task with periodical speaker-listener changes,
and we would like to evaluate ISC in both directions from segments with the same speaker-listener roles.

In this script we generate random coupled data ("speaker" and "listener") with known coefficients and SNR level
and evaluate ISC with different segmentation strategies.

The following options are considered:
- Baseline: Performance of ISC for a long, continuous recordings with fixed coupling coefficients (beta).
- Concatenation: Selected segments are concatenated and then ISC is estimated from the concatenated data.
- Averaged-over-segments: ISC is estimated separately for segments and the results are then averaged.

A summary plot is generated at the end.
'''


import numpy as np
import sys
import matplotlib as mplt
import matplotlib.pyplot as plt
# add module from "hypercodes" git repo
sys.path.append('/home/adamb/hypercodes/current_code/')
import couplingISC as isc


# Parameters for random data generation
tr_no = 1500  # keep it a multiple of 30
vox_no = 1000
data_dims = (tr_no, vox_no)
shift = 3  # no. of TRs to consider in both directions when calculating coupling ISC
snr_levels = np.logspace(-1, 1, 40)  # 40 log-spaced values from 0.1 to 10
beta = np.asarray([0, 0.2, 0, 0.2, 0.5, 2, 0])  # coefficients, shift * 2 + 1 values
beta = beta/np.sum(beta)  # normalize coefficients so that the generated coupled data has the same scale

# Preallocate vars holding ISC results for all SNR levels.
# Column 1 for holding the mean error of coefficients
# Column 2 for holding the STD of error across voxels
baseline_res = np.zeros((len(snr_levels), 2))  # baseline ISC on long continuous coupled data
concat_res = np.zeros((len(snr_levels), 2))  # concatenated independent segments
aver_segments_res = np.zeros((len(snr_levels), 2))  # ISC on segments, then averaging over segment-level results

# params for segmentation
bounds = np.arange(0, tr_no, tr_no/30)  # thirty segments, every second will be concatenated
bounds1 = bounds[0::2].astype(int)  # every second segment will be selected for concatenation
bounds2 = bounds[1::2].astype(int)

print('\n\nCoupling ISC estimation.')

# iterate over snr levels
for snr_it, snr_level in enumerate(snr_levels):

    # progress tracking
    if snr_it % 5 == 0:
        print('SNR level: ' + str(snr_level))

    # get noise level from SNR level - the former is needed for get_coupled_set_2d()
    noise_level = 1/snr_level

    # generate random data
    X, Y = isc.get_coupled_set_2d(data_dims, beta, shift, noise_level)

    #################################################################
    # Baseline ISC on long continuous data
    #################################################################
    # calculate baseline ISC
    beta_est = isc.coupling_onestep(X[1:tr_no//2, :], Y[1:tr_no//2, :], shift)  # on only half of the data

    # get difference in betas (~ accuracy)
    baseline_res[snr_it, 0] = np.mean(np.sqrt(np.sum((beta_est - beta) ** 2, 1)))  # mean sqrt of sum of squared diffs
    baseline_res[snr_it, 1] = np.std(np.sqrt(np.sum((beta_est - beta) ** 2, 1))) # std of the same

    #################################################################
    # ISC on concatenated independent segments
    #################################################################
    # concatenate segments
    for i in range(len(bounds1)):
        if i == 0:
            X1 = X[bounds1[i]:bounds2[i], :]
            Y1 = Y[bounds1[i]:bounds2[i], :]
        else:
            X1 = np.concatenate((X1, X[bounds1[i]:bounds2[i], :]), axis=0)
            Y1 = np.concatenate((Y1, Y[bounds1[i]:bounds2[i], :]), axis=0)

    # calculate ISC
    beta_est = isc.coupling_onestep(X1, Y1, shift)

    # get difference in betas (~ accuracy)
    concat_res[snr_it, 0] = np.mean(np.sqrt(np.sum((beta_est - beta) ** 2, 1)))  # mean sqrt of sum of squared diffs
    concat_res[snr_it, 1] = np.std(np.sqrt(np.sum((beta_est - beta) ** 2, 1)))  # std of the same

    #################################################################
    # ISC on independent segments followed by averaging
    #################################################################
    # preallocate var holding segment-level betas
    segments_beta_est = np.zeros((len(bounds1), vox_no, shift*2+1))

    # iterate over segments
    for i in range(len(bounds1)):
        X1 = X[bounds1[i]:bounds2[i], :]
        Y1 = Y[bounds1[i]:bounds2[i], :]

        # calculate ISC
        segments_beta_est[i, :, :] = isc.coupling_onestep(X1, Y1, shift)

    # average over the results of segment-level ISCs
    beta_est = np.mean(segments_beta_est, 0)

    # get difference in betas (~ accuracy)
    aver_segments_res[snr_it, 0] = np.mean(np.sqrt(np.sum((beta_est - beta) ** 2, 1)))  # mean sqrt of sum of squared diffs
    aver_segments_res[snr_it, 1] = np.std(np.sqrt(np.sum((beta_est - beta) ** 2, 1)))  # std of the same


# plot the results
errbarCont = plt.errorbar(snr_levels, baseline_res[:, 0], yerr=baseline_res[:, 1], fmt='bo-')
plt.errorbar(snr_levels, concat_res[:, 0], yerr=concat_res[:, 1], fmt='ro-')
plt.errorbar(snr_levels, aver_segments_res[:, 0], yerr=aver_segments_res[:, 1], fmt='go-')
errbarCont[0].axes.set_xscale('log')
errbarCont[0].axes.set_xticks([0.1, 1, 10])
errbarCont[0].axes.get_xaxis().set_major_formatter(mplt.ticker.ScalarFormatter())
plt.title('Coeff detection error as a function of SNR\n', fontsize=20)
plt.xlabel('\nSNR', fontsize=18)
plt.ylabel('Mean summed error of coeffs\n', fontsize=18)
plt.legend(labels=['Baseline (ground truth)', 'Concatenated segments', 'Averaged over segments'])
plt.rcParams.update({'font.size': 18})













