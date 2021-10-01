#!/usr/bin/env conda run -n py27Env python2.7
# -*- coding: utf-8 -*-

"""
Hyperalign on one half of a hyperscanning task and look
for improvements in leave-one-out ISC in the other half.
"""

import numpy as np
import pickle
from mvpa2.suite import *
import time
import glob
import sys
sys.path.append('/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/support_scripts/')
from phaseScramble_2 import *
from scipy import stats

def main():

    print('\nlets hyperalign\n')

    # define hyperscanning task descriptions
    taskDescrips = ['storytelling_independent',
                    'storytelling_joint',
                    'listening',
                    'reading']

    # parameters
    debug = False
    task = 4 # see task descriptions above
    radius = 3 # number of voxels in hyperalignment searchlight radius
    sparse_radius = 3 # number of voxels between neighboring searchlight spheres
    nproc = 10 # number of parallel processes to feed into hyperalignment function

    # set dataset labels
    dsLabels = ['train','test']

    # set base folder
    baseFolder = '/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/preprocessing/hyperalignment/'

    if debug:
        datasetFile = baseFolder + 'datasets/debug_' + taskDescrips[task]
    else:
        datasetFile = baseFolder + 'datasets/' + taskDescrips[task]

    # load training and testing data
    ds_all = h5load(datasetFile)

    # get training and testing sample indices (half and half)
    order = 0 # 0 = train on first half, test on second, 1 = the opposite
    halfSampleNum = np.round(ds_all[0].nsamples / 2)
    sampleInds = [[]] * 2
    if order == 0:
        sampleInds[0] = np.arange(halfSampleNum) # training sample indices
        sampleInds[1] = np.arange(halfSampleNum,ds_all[0].nsamples,1) # testing sample indices
    else:
        sampleInds[0] = np.arange(halfSampleNum, ds_all[0].nsamples, 1) # training sample indices
        sampleInds[1] = np.arange(halfSampleNum) # testing sample indices

    # get number of subjects in full dataset
    numSubs = len(ds_all)

    # split up into training and testing datasets
    ds = [[]] * 2 # initialize
    for DS in range(len(ds)): # for each data set (0=training, 1=testing)
        ds[DS] = [[]] * numSubs # initialize
        for SUB in range(numSubs): # for each subject
            ds[DS][SUB] = ds_all[SUB][sampleInds[DS],:]
            ds[DS][SUB].samples = stats.zscore(ds[DS][SUB].samples, axis=0)

    # verify that subject ID lists are identical between training and testing sets

    # for each dataset...
    EPIdata = [[]] * 2
    corrData = [[]] * 2
    medCorr = [[]] * 2
    for DS in range(2):

        # get number of subjects
        numSubs = len(ds[DS])

        # get EPI dimensions (samples x voxels)
        dims = np.array(ds[DS][0].__array__().shape)

        # initialize raw EPI data array
        EPIdata[DS] = np.empty([dims[0], dims[1], len(ds[DS])])

        # initialize pre-hyperalignment ISC coefficient array (subs x voxels)
        corrData[DS] = np.empty([numSubs, dims[1]])

        # for each subject...
        for SUB in range(numSubs):

            # get EPI data
            EPIdata[DS][:,:,SUB] = ds[DS][SUB].__array__()

        # for each subject...
        for SUB in range(numSubs):

            # get mean of data from all participants EXCEPT the current participant
            otherSubs = np.arange(0, numSubs)
            otherSubs = np.delete(otherSubs, SUB)
            groupMean = np.mean(EPIdata[DS][:,:,otherSubs], axis=2)

            # get correlation between current participant and groupMean
            corrData[DS][SUB, :] = fastColumnCorr(EPIdata[DS][:,:,SUB], groupMean)

        # get median ISC across participants
        medCorr[DS] = np.median(corrData[DS], axis=0)
        print('mean (across voxels) median (across subs) corr in ' + dsLabels[DS] + ' set BEFORE hyperalignment: ' + str(np.round(np.mean(medCorr[DS]),3)))

    # we call SearchlightHyperalignment mostly with default values:
    # each sphere has a radius of 3 voxels, sphere centers are also 3 voxels apart,
    # all voxels in a given sphere are used for alignment
    slhyper = SearchlightHyperalignment(radius=radius,
                                        sparse_radius=sparse_radius,
                                        nproc=nproc)

    # call the hyperalignment object with the full dataset we have,
    # resulting mappers will be stored in slhypmaps
    slhyperStart = time.time()
    slhypmaps = slhyper(ds[0])
    print('\nHyperalignment took ' + str(time.time() - slhyperStart) + ' secs')

    # compute post-hyperalignment metrics
    ds_hyper = [[]] * 2
    EPIdata_hyper = [[]] * 2
    corrData_hyper = [[]] * 2
    medCorr_hyper = [[]] * 2
    for DS in range(2):

        # Applying hyperalignment parameters is similar to applying any mapper in
        # PyMVPA. We apply the hyperalignment parameters by running the dataset
        # through the forward() function of the mapper.
        ds_hyper[DS] = [h.forward(sd) for h, sd in zip(slhypmaps, ds[DS])]

        # get EPI dimensions (samples x voxels)
        dims = np.array(ds_hyper[DS][0].__array__().shape)

        # initialize raw EPI data array
        EPIdata_hyper[DS] = np.empty([dims[0], dims[1], len(ds_hyper[DS])])

        # initialize pre-hyperalignment ISC coefficient array (subs x voxels)
        corrData_hyper[DS] = np.empty([numSubs, dims[1]])

        # for each subject...
        for SUB in range(numSubs):
            # get EPI data
            EPIdata_hyper[DS][:, :, SUB] = ds_hyper[DS][SUB].__array__()

        # for each subject...
        for SUB in range(numSubs):
            # get mean of data from all participants EXCEPT the current participant
            otherSubs = np.arange(0, numSubs)
            otherSubs = np.delete(otherSubs, SUB)
            groupMean = np.mean(EPIdata_hyper[DS][:, :, otherSubs], axis=2)

            # get correlation between current participant and groupMean
            corrData_hyper[DS][SUB, :] = fastColumnCorr(EPIdata_hyper[DS][:, :, SUB], groupMean)

        # get median ISC across participants
        medCorr_hyper[DS] = np.median(corrData_hyper[DS], axis=0)
        print('mean (across voxels) median (across subs) corr in ' + dsLabels[
            DS] + ' set BEFORE hyperalignment: ' + str(np.round(np.mean(medCorr_hyper[DS]), 3)))

    # save name
    if task == 3:
        saveFile = baseFolder + 'results/listening_5050'
    else:
        saveFile = baseFolder + 'results/reading_5050'
    if debug:
        saveFile = saveFile + '_debug'
    print('saving files to: ')
    print(saveFile + '_med_corr_pre_hyp')
    print(saveFile + '_med_corr_post_hyp')

    # save median correlation array
    np.save(saveFile + '_med_corr_pre_hyp', medCorr)
    np.save(saveFile + '_med_corr_post_hyp',medCorr_hyper)
    # print('saving output...')
    # with open(saveFile + '.pkl', 'wb') as f:
    #     pickle.dump([medCorr, medCorr_hyper], f, protocol=2)
    # print('output saved here: ' + saveFile + '.pkl')

    # save mapping
    h5save(saveFile + '_hyperMappings.hd5', slhypmaps)
    print('yay')

if __name__ == '__main__':
    main()
