#!/usr/bin/env conda run -n py27Env python2.7
# -*- coding: utf-8 -*-

"""
Hyperalign on control dataset (either listening or reading or both)
and apply transformation to storytelling task and look for improvements
in leave-one-out-ISC.
Input parameters:
    --trainTasktrain (-r): control task data on which to hyperalign
        1=storytelling_independent
        2=storytelling_joint
        3=listening
        4=reading
        5=listening_&_reading
    --testCond (-c): storytelling condition on which to test
        0=independent
        1=joint
    --testTask (-e): subset (if any) of storytelling task on which to test
        0=all (no subsetting)
        1=storytelling listening intervals
        2=storytelling reading intervals
"""

import numpy as np
import pickle
from mvpa2.suite import *
from scipy import stats
import time
import glob
import sys
sys.path.append('/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/support_scripts/')
from phaseScramble_2 import *

def main():

    # get pumped
    print('\nlets hyperalign\n')

    #####################################
    ### housekeeping ####################
    #####################################

    # initialize parser
    parser = argparse.ArgumentParser()

    # parse arguments...
    parser.add_argument(
        '-r',
        '--trainTask',
        default=5,
        type=int,
        help='1=storytelling_independent, 2=storytelling_joint, 3=listening, 4=reading, 5=listening_&_reading')

    parser.add_argument(
        '-c',
        '--testCond',
        default=0,
        type=int,
        help='0=independent, 1=joint')

    parser.add_argument(
        '-e',
        '--testTask',
        default=0,
        type=int,
        help='0=all, 1=storytelling listening intervals, 2=storytelling reading intervals')

    # parse arguments, get list
    args = parser.parse_args()

    #####################################
    ### housekeeping ####################
    #####################################

    # define hyperscanning task descriptions
    taskDescrips = ['storytelling_independent',
                    'storytelling_joint',
                    'listening',
                    'reading',
                    'listening_&_reading']

    # parameters
    trainTask = args.trainTask # see descriptions above *use 1-5 here, not 0-4
    testCond = args.testCond # 0=independent, 1=joint
    testTask = args.testTask # 0=all, 1=storytelling listening intervals, 2=storytelling reading intervals
    radius = 3  # number of voxels in hyperalignment searchlight radius
    sparse_radius = 3  # number of voxels between neighboring searchlight spheres
    nproc = 10  # number of parallel processes to feed into hyperalignment function

    # save file
    baseFolder = '/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/preprocessing/hyperalignment/'
    saveFile = baseFolder + 'results/train_on_' + taskDescrips[trainTask-1] + '_test_on_' + taskDescrips[testCond]
    if testTask == 0:
        testTag = 'all'
    elif testTask == 1:
        testTag = 'listening'
    elif testTask == 2:
        testTag = 'speaking'
    saveFile = saveFile + '_' + testTag

    # feedback
    print('training on ' + taskDescrips[trainTask - 1] + ', testing on ' + taskDescrips[testCond] + ' ' + testTag)

    # set dataset labels
    dsLabels = ['train','test']

    #####################################
    ### prepare datasets ################
    #####################################

    # get dataset paths (0=train, 1=test)
    if trainTask == 5:
        trainFiles = [baseFolder + 'datasets/listening',
                      baseFolder + 'datasets/reading']
    else:
        trainFiles = [baseFolder + 'datasets/' + taskDescrips[trainTask-1]]
    testFile = [baseFolder + 'datasets/' + taskDescrips[testCond]]
    datasetFiles = trainFiles + testFile

    # load source datasets
    dsSource = [[]] * len(datasetFiles)
    for DS in range(len(dsSource)):
        dsSource[DS] = h5load(datasetFiles[DS])

    # get reference subject list
    refList = dsSource[0][0].a['full_dataset_subID_list'].value

    # verify that subject ID lists are identical between training and testing sets
    # if not identical, rearrange datasets to line them up, arbitrarily using first
    # dataset as reference (see above)
    for DS in np.arange(1,len(dsSource)): # for datasets 2 on...

        # get current dataset's subject list
        subList = dsSource[DS][0].a['full_dataset_subID_list'].value

        # if it's different from the reference list...
        if subList != refList:

            if len(subList) != len(refList):
                print('Warning! Number of subjects does not match between datasets!')
                break
            else:

                # preallocate reordering indices
                reorderInds = np.empty(len(subList),dtype=int)

                # for each subject on the reference list...
                for SUB in range(len(refList)):

                    # get the corresponding index in the current subject list
                    reorderInds[SUB] = subList.index(refList[SUB])

                # make a temporary copy of the current dataset
                tmp = dsSource[DS][:]

                # for each subject on the reference list...
                for SUB in range(len(refList)):

                    # add the reordered the subject list to the current subject's dataset
                    dsSource[DS][SUB].a['full_dataset_subID_list_REORDERED'] = [subList[SUB2] for SUB2 in reorderInds]

                    # update the dataset at the current position per reorderInds
                    dsSource[DS][SUB] = tmp[reorderInds[SUB]]

                    # tmp check
                    print('\nchecking reordered datasets...')
                    print('ref sub ' + str(SUB + 1) + ' ID: ' + dsSource[0][SUB].a['subID'].value)
                    print('cur sub ' + str(SUB + 1) + ' ID: ' + dsSource[DS][SUB].a['subID'].value)
                    if dsSource[0][SUB].a['subID'].value != dsSource[DS][SUB].a['subID'].value:
                        print('something went wrong!')

    # preallocate list for training and testing data
    ds = [[]] * 2
    for DS in [0,1]:
        ds[DS] = [[]] * len(dsSource[0])

    # if we're training on both the listening and reading tasks...
    if trainTask == 5:

        # for each subject in the listening data set...
        for SUB in range(len(dsSource[0])):

            # vertically stack listening and reading datasets and add to ds
            ds[0][SUB] = vstack([dsSource[0][SUB],dsSource[1][SUB]])
    else:
        ds[0] = dsSource[0]

    # add test data to ds
    ds[1] = dsSource[len(dsSource)-1]

    # set number of samples that should be in each storytelling dataset
    totTRs = 1230

    # check that the storytelling datasets have the proper number of samples
    if ds[1][0].samples.shape[0] != totTRs:
        print('WARNING! Storytelling data does not have the expected number of samples!')

    # number of TRs per turn
    TRsPerTurn = 41

    # number of speech turns per participant per condition
    numTurns = int(totTRs / TRsPerTurn)

    # get speaker/listener TR indices
    turnTRs = [[]] * numTurns
    for TURN in range(int(numTurns)):
        if TURN == 0:
            inds = np.array(list(range(TRsPerTurn)))
        else:
            inds = inds + TRsPerTurn
        turnTRs[TURN] = inds.reshape([len(inds),1])

    # get speaker masks
    trSpeakerMask = [[]] * 2
    trSlicer = [[]] * 2
    for SITE in [0, 1]:  # 0=cbs,1=dbic
        seed = SITE
        trSpeakerMask[SITE] = [[]] * int(numTurns / 2)
        for TURN in range(int(numTurns/2)):
            trSpeakerMask[SITE][TURN] = seed
            seed += 2
        trSlicer[SITE] = [turnTRs[i] for i in trSpeakerMask[SITE]]

    # Conditionally slice test data to only include the relevant speech turn type
    # (listening or speaking). CBS participants spoke first in all cases.
    if testTask > 0:

        # for each subject in the storytelling data...
        for SUB in range(len(ds[1])):

            # if it's a cbs participant...
            if ds[1][SUB].a['subID'].value[0] == 'h':

                # if we're testing on listener turns
                if testTask == 1: # if testing on listening task

                    # get listening TRs
                    testTRs = np.vstack(trSlicer[1]).reshape(-1)

                # if testing on speaker turns...
                elif trainTask == 4:

                    # get reading TRs
                    testTRs = np.vstack(trSlicer[0]).reshape(-1)

            # if it's a dbic participant...
            elif ds[1][SUB].a['subID'].value[0] == 's':

                # if we're testing on listener turns
                if testTask == 1:  # if testing on listening task

                    # get listening TRs
                    testTRs = np.vstack(trSlicer[0]).reshape(-1)

                # if testing on speaker turns...
                elif trainTask == 4:

                    # get reading TRs
                    testTRs = np.vstack(trSlicer[1]).reshape(-1)

            # slice the dataset
            ds[1][SUB] = ds[1][SUB][testTRs,:]

    # standardize data
    for DS in [0, 1]:
        for SUB in range(len(ds[DS])):
            ds[DS][SUB].samples = stats.zscore(ds[DS][SUB].samples, axis=0)

    # tmp CHECK
    print(str(ds[0][0].samples.shape[0]) + ' samples in test dataset')

    ######################################################
    ### compute ISC before hyperalignment ################
    ######################################################

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

    ########################################################
    ### HYPERALIGN! ########################################
    ########################################################

    # we call SearchlightHyperalignment mostly with default values:
    # each sphere has a radius of 3 voxels, sphere centers are also 3 voxels apart,
    # all voxels in a given sphere are used for alignment
    slhyper = SearchlightHyperalignment(radius=radius,
                                        sparse_radius=sparse_radius,
                                        nproc=nproc)

    # call the hyperalignment object with the full dataset we have,
    # resulting mappers will be stored in slhypmaps
    print('\ncommencing hyperalignment...')
    slhyperStart = time.time()
    slhypmaps = slhyper(ds[0])
    print('\nHyperalignment took ' + str(time.time() - slhyperStart) + ' secs')

    #########################################
    ### recompute ISC in the common space ###
    #########################################

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

    # save median correlation array
    np.save(saveFile + '_med_corr_pre_hyp', medCorr)
    np.save(saveFile + '_med_corr_post_hyp',medCorr_hyper)

    # save mapping
    h5save(saveFile + '_hyperMappings.hd5', slhypmaps)

    print('\n\n\nThe loony bin is closing shop, greetings',
          'from The Combine!\n\n')

if __name__ == '__main__':
    main()
