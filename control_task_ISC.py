#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script runs intersubject correlation (ISC) analyses on the hyperscanning
control task (listening and reading) data. It first loads all of the relevant
time series, then computes the correlation between each individual's time
series and the mean time series across the rest of the group. It then generates
a null distribution for each voxel in each subject by iteratively scrambling the
subject's time series (either through phase scrambling or circle shifting) and
recomputing the correlation with the rest of the group. The true correlation
value is then evaluated against the null distribution to get a p-value.
P-values are FDR corrected.
Major parameters of note:
    debug: boolean, if True, truncates all time series to a small, manageable size (debugTRs x debugVox) to speed things up
    permutations: number of iterations in permutation test
    parallel: boolean, if True, use the joblib parallelization function (speeds things up considerably)
    numJobs: number of parallel processes for joblib to implement
    verbosity: joblib verbosity -- don't go over 50 lest ye want to print like a million outputs and slow everything down
    circShift: boolean, if True, use the circle shift method to "scramble" time series, if False, use phase scrambling
    fisherZtran: boolean, if True, fisher z transform correlation coefficients
    normalize: boolean, if True, z-score time series before computing correlations
    firProc: numeric that indicates which fitting function to use to extract null distribution parameters (see below)
    alpha: alpha to use for permutation tests
    twoTailed: boolean, if True, uses a two-tailed permutation test, if False, uses a right-tailed test
    detrend: boolean, if True, uses epi timeseries that were preprocessed with the nilearn detrending input set to true
    saveOutput: boolean, if True, saves some ISC output variables in a .pkl file
"""

# imports
import argparse
import scipy.io as sio
import numpy as np
import pandas as pd
import time
from joblib import Parallel, delayed
from scipy import stats
import statsmodels.stats.multitest as multi
from sklearn import preprocessing
import pickle
from distutils import util
import sys
import math
sys.path.append('/dartfs-hpc/rc/home/z/f00589z/hyperscanning/support_scripts/')
from phaseScramble import *
from CircleShift import *


# parallel subject loop wrapper function
def parallelSubWrapper(subNum,numSubs,boldData,permutations,circShift,fisherZtran,twoTailed,alpha):
    """
    function to use with joblib below to run permutation tests for participants in parallel
    :param subNum: participant index
    :param numSubs: total number of participants in the analysis
    :param boldData: timeseries data [timepoints x voxels]
    :param permutations: number of permutations to use
    :param circShift: boolean indicating whether or not (True or False, respectively) to use the circle shift method to "scramble" timeseries
    :param fisherZtran: boolean indicating whether or not fisher z transform correlation coefficients
    :param twoTailed: boolean indicating whether to use a two tailed (True) or one tailed (False) permutation test
    :return:
        corrData: real voxelwise correlations between participant and the mean across all other participants [voxels x 1]
        nullCorr: null voxelwise correlations between participant and the mean across all other participants [permutations x voxels]
        permTest: array containing various results from permutation test and normal fits to the null distribution (see permTestMap variable above and comments below)
            structure of permTest:
            permTest[0]: permutation test p-values [voxels x 1]
            permTest[1]: FDR corrected permutation test p-values [voxels x 1]
            permTest[2]: p-value summary stats
                permTest[2][0]: number of voxels with uncorrected p-vals < alpha
                permTest[2][1]: proportion of voxels with uncorrected p-vals < alpha
                permTest[2][2]: number of voxels with FDR corrected p-vals < alpha
                permTest[2][3]: proportion of voxels with FDR corrected p-vals < alpha
            permTest[3]: null dist normal fit parameters [voxels x 1]
                permTest[3][voxel][0]: mean
                permTest[3][voxel][1]: standard deviation
            permTest[4]: Kolmogorov–Smirnov goodness of fit test results [voxels x 1]
                permTest[4][voxel][0]: KS test statistic
                permTest[4][voxel][1]: p-value (values less than alpha suggest poor fit)
            permTest[5]: logical vector for KS goodness of fit where 0 = good fit, 1 = bad fit [voxels x 1]
            permTest[6]: table indicating number and proportion of bad normal fits to the perm-based null distribution [1 x 2 table]
    """

    # map of what each participant-specific sublist in permTestMap contains
    permTestMap = ['permPval', 'FDR_corrected_permPval', 'propSigFDRvoxels', 'normParams', 'KStest', 'badFits',
                   'badFitsSummary']

    # get mean of data from all participants EXCEPT the current participant
    otherSubs = np.arange(0,numSubs)
    otherSubs = np.delete(otherSubs,subNum)
    ogMeanMethod = False;
    if ogMeanMethod:

        # original method for computing the mean across remaining subjects for ISC
        groupMean = np.mean([boldData[i] for i in otherSubs], axis=0)

    else: # new way (to avoid numpy memory errors)

        # new method attempting to the memory burden on numpy -- breaks the mean computation down into chunks
        chunkLength = 100 # set number of chunks (we will divide the number of voxels into this many chunks)
        numVox = boldData[0].shape[1] # get number of voxels in bold data
        numChunks = math.ceil(numVox / chunkLength) # get chunk length -- number of voxels in current chunk
        startVox = 0  # initialize starting voxel per chunk counter
        groupMean = np.empty(boldData[0].shape)
        breakFlag = False
        for CHUNK in range(numChunks):

            endVox = startVox + chunkLength
            if endVox > numVox:
                endVox = numVox
                breakFlag = True
            chunkInds = np.arange(startVox, endVox)
            groupMean[:,chunkInds] = np.mean([boldData[i][:,chunkInds] for i in otherSubs], axis=0)

            # get starting voxel for next chunk
            startVox += chunkLength

            if breakFlag:
                break

    # get REAL correlation between current participant and groupMean
    corrData = fastColumnCorr(boldData[subNum], groupMean)

    # fisher z transform if selected
    if fisherZtran:
        corrData = np.arctanh(corrData)

    # generate null correlations
    nullCorr = [[]] * permutations
    feedbackMultiple = round(permutations / 10) # get number to use as a flag to generate feedback roughly 10 times while null correlations are being generated
    for PERM in range(permutations): # for each permutation...
        if circShift:
            nullCorr[PERM] = fastColumnCorr(Circle_Shift(boldData[subNum]),groupMean)
        else:
            nullCorr[PERM] = fastColumnCorr(phase_scrambling(boldData[subNum]),groupMean)
        if PERM % feedbackMultiple == 0: # if the current permutation is zero or a multiple of feedbackMultiple
            print('getting null correlation for subject ' + str(subNum + 1) + ', permutation ' + str(PERM + 1)) #provide feedback about which subject and which permutation we're on

    # convert permutation data to a permutations x voxels array
    nullCorr = np.asarray(nullCorr)

    # fisher z transform if selected
    if fisherZtran:
        nullCorr = np.arctanh(nullCorr)

    # preallocate sublists of permTest array
    permTest = [[]] * len(permTestMap)

    # get permutation test p-values and apply FDR correction
    permTest[0] = [[]] * nullCorr.shape[1] # permutation test p-values
    for VOX in range(nullCorr.shape[1]): # for each voxel...
        if twoTailed:
            permTest[0][VOX] = len(np.where(abs(nullCorr[:,VOX]) > abs(corrData[VOX]))[0]) / float(permutations) # proportion of permutations where absolute value of null correlation is greater than absolute value of real correlation
            alphaPrime = alpha / 2
        else:
            permTest[0][VOX] = len(np.where(nullCorr[:,VOX] > corrData[VOX])[0]) / float(permutations) # proportion of permutations where value of null correlation is greater than value of real correlation
            alphaPrime = alpha
    permTest[1] = multi.fdrcorrection(permTest[0], alpha = alphaPrime) # FDR correction - logical vector with length = # voxels
    permTest[2] = [[]] * 4 # preallocate lists for permutation test summary info
    permTest[2][0] = len(np.where(np.array(permTest[0]) < alphaPrime)[0]) # number of voxels that show significant correlations with the group
    permTest[2][1] = permTest[2][0] / nullCorr.shape[1] # proportion of voxels that show significant correlations with the group
    permTest[2][2] = np.count_nonzero(permTest[1][0]) # number of voxels that show significant FDR CORRECTED correlations with the group
    permTest[2][3] = len(np.where(permTest[1][0])[0]) / nullCorr.shape[1] # proportion of voxels that show significant FDR CORRECTED correlations with the group

    # fit normal distributions to null data
    permTest[3] = [[]] * nullCorr.shape[1] # preallocate list for voxelwise normal dist parameters
    permTest[4] = [[]] * nullCorr.shape[1] # preallocate list for voxelwise Kolmogorov–Smirnov test results
    permTest[5] = np.zeros((nullCorr.shape[1],), dtype=int) # preallocate voxelwise array to indicate bad fits
    for VOX in range(nullCorr.shape[1]): # for each voxel...

        # fit normal distribution - using 3 different methods here as a sanity check
        permTest[3][VOX] = [[]] * 2 # preallocate for mean and SD
        permTest[3][VOX][0] = np.mean(nullCorr[:,VOX])
        permTest[3][VOX][1] = np.std(nullCorr[:,VOX])
        permTest[4][VOX] = stats.kstest(nullCorr[:,VOX], "norm", permTest[3][VOX]) # measure goodness of fit
        if permTest[4][VOX][1] < alpha: # if a voxel has a bad normal fit...
            permTest[5][VOX] = 1 # flag it
    permTest[6] = pd.DataFrame(data={'numBadFits': [sum(permTest[5])], 'propBadFits': [sum(permTest[5]) / len(permTest[5])]}) # number and proportion of voxels with bad normal fits to the null distribution

    return corrData, nullCorr, permTest

#%% main
def main():

    """
    Note that a lot of the inputs to argparse should really be
    booleans, but because thise script is designed to take inputs
    from the command line, True and False inputs are instead
    inputted as strings, then converted to booleans later.
    :return:
    """
    # parse ipnut arguments
    parser = argparse.ArgumentParser()

    # Input arg "folder"
    parser.add_argument(
        '-n',
        '--numNodes',
        default=16,
        type=int,
        help='number of requested nodes')

    parser.add_argument(
        '-c',
        '--circShift',
        default='True',
        type=str,
        help='True = use circle shifting, False = use phase scrambling')

    parser.add_argument(
        '-d',
        '--debug',
        default='False',
        type=str,
        help='True = run in debug mode')

    parser.add_argument(
        '-p',
        '--permutations',
        default=1000,
        type=int,
        help='number of permutations')

    parser.add_argument(
        '-l',
        '--parallel',
        default='True',
        type=str,
        help='True = use joblib to run permutation tests in parallel')

    parser.add_argument(
        '-v',
        '--verbosity',
        default=50,
        type=int,
        help='verbosity of joblib -- be warned that going above 50 results in a wild amount of output'
    )

    parser.add_argument(
        '-f',
        '--fisherZtran',
        default='True',
        type=str,
        help='True = fishers z transform null distributions')

    parser.add_argument(
        '-o',
        '--normalize',
        default='True',
        type=str,
        help='z score EPI time series before ISC')

    parser.add_argument(
        '-a',
        '--alpha',
        default=0.05,
        type=float,
        help='alpha for permutation tests')

    parser.add_argument(
        '-w',
        '--twoTailed',
        default='True',
        type=str,
        help='True = two tailed permutation tests, False = one tailed')

    parser.add_argument(
        '-e',
        '--detrend',
        default='False',
        type=str,
        help='True = use detrended EPI time series')

    parser.add_argument(
        '-s',
        '--saveOutput',
        default='True',
        type=str,
        help='True = save ISC results')

    # parse arguments, get list
    args = parser.parse_args()

    # paths
    baseFolder = '/dartfs-hpc/rc/home/z/f00589z/hyperscanning/'
    inputFolder = baseFolder + 'control_tasks/nuisRegr_output_files/'
    outputFolder = baseFolder + 'control_tasks/control_ISC_output/'
    subListFolder = baseFolder + 'misc/'

    # mark starting time
    startTime = time.time()

    # load hyperscanning subject list
    subList = pd.read_pickle(subListFolder + 'hyperscanning_subject_list.pkl')

    # get number of participants
    numSubs = subList.shape[0]

    # define task names
    taskNames = np.array(['listening','reading'])

    # initialize the saveTag variable
    saveName = 'controlISC'

    # set major parameters

    # indicate whether or not we're debugging (if so, use a small subset of TRs and voxels to speed things up)
    debug = bool(util.strtobool(args.debug))
    debugTRs = 201
    debugVox = 800
    if debug:
        saveName = saveName + '_debug'

    # set number of permutations
    permutations = args.permutations
    saveName = saveName + '_' + str(permutations) + 'perm'

    # use joblib or not
    parallel = bool(util.strtobool(args.parallel))
    if not parallel:
        saveName = saveName + '_notParallel'

    # number of jobs for joblib
    numJobs = args.numNodes

    # set joblib verbosity -- don't go over 50 lest ye want to print like a million outputs and slow everything down
    verbosity = args.verbosity

    # select whether or not to use the circle shifting method
    circShift = bool(util.strtobool(args.circShift))
    if circShift:
        saveName = saveName + '_cShift'
    else:
        saveName = saveName + '_pScram'

    # select whether or not to Fisher transform the resulting pearson correlation distributions
    fisherZtran = bool(util.strtobool(args.fisherZtran))
    if fisherZtran:
        saveName = saveName + '_fishZ'

    # select whether or not to z-scale timeseries before ISC
    normalize = bool(util.strtobool(args.normalize))
    if normalize:
        saveName = saveName + '_norm'

    # set alpha for permutation tests
    alpha = args.alpha

    # two or one-tailed (right tail) permutation test
    twoTailed = bool(util.strtobool(args.twoTailed))
    if twoTailed:
        saveName = saveName + '_twoTailed'
    else:
        saveName = saveName + '_oneTailed'

    # use detrended data?
    useDetrendedData = bool(util.strtobool(args.detrend))
    if useDetrendedData:
        loadTag = 'detrended'
        saveName = saveName + '_detrended'

    # select whether or not to save output
    saveOutput = bool(util.strtobool(args.saveOutput))

    ##############################################
    ### load listening and reading time series ###
    ##############################################

    # loop through participants...
    boldData = [[]] * 2
    for TASK in [0,1]: # for each task, listening, then reading...

        # preallocate task data list
        boldData[TASK] = [[]] * numSubs

        for SUB in range(numSubs):

            # get file name
            fileName = inputFolder + 'sub-' + subList['subID'][SUB] + '_ses-pair0' + str(subList['pairNum'][SUB]) + '_task-storytelling' + str(TASK + 3) + '_run-0' + str(TASK + 3) + '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021_interp.mat'

            # load data
            tmp = sio.loadmat(fileName) # load file
            boldData[TASK][SUB] = tmp['tseries'] # get time series data
            print('\nloaded ' + str(boldData[TASK][SUB].shape[0]) + ' x ' + str(boldData[TASK][SUB].shape[1]) + ' timeseries for ' + taskNames[TASK] + ' task, sub ' + subList['subID'][SUB])

            if normalize:
                boldData[TASK][SUB] = preprocessing.normalize(boldData[TASK][SUB])
                print('normalizing timeseries')

            if debug:

                # take a small subset of each time series to speed up the debugging process
                boldData[TASK][SUB] = boldData[TASK][SUB][np.ix_(np.arange(debugTRs),range(debugVox))]

    #########################
    ### permutation tests ###
    #########################

    # preallocate correlation lists
    corrData = [[]] * 2
    nullCorr = [[]] * 2
    permTest = [[]] * 2

    # Voxelwise correlation between participant and the rest of the group (mean)
    for TASK in [0,1]: # for each task...

        corrData[TASK] = [[]] * numSubs
        nullCorr[TASK] = [[]] * numSubs
        permTest[TASK] = [[]] * numSubs

        # run ISC
        if parallel:

            # run joblib
            tmp = Parallel(n_jobs=numJobs, verbose=verbosity)(delayed(parallelSubWrapper)(SUB,numSubs,boldData[TASK],permutations,circShift,fisherZtran,twoTailed,alpha) for SUB in range(numSubs))

            # assign joblib outputs
            for SUB in range(numSubs):
                corrData[TASK][SUB] = tmp[SUB][0]
                nullCorr[TASK][SUB] = tmp[SUB][1]
                permTest[TASK][SUB] = tmp[SUB][2]
        else:

            for SUB in range(numSubs):
                corrData[TASK][SUB], nullCorr[TASK][SUB], permTest[TASK][SUB] = parallelSubWrapper(SUB,numSubs,boldData[TASK],permutations,circShift,fisherZtran,twoTailed,alpha)

        # print some permutation test and goodness of fit summary info
        print('\nFinished processing participant ' + str(SUB + 1) + ' of ' + str(numSubs) + ' for task ' + str(TASK + 3))
        print('% voxels with FDR corrected significant correlation with group: ' + str((permTest[TASK][SUB][2][1]) * 100) + '%')
        print('% voxels for which null dist. was normal: ' + str((1 - permTest[TASK][SUB][6].iloc[0,1]) * 100) + '%')

    #########################################################
    ### get group level null distributions for each voxel ###
    #########################################################

    if numSubs > 1: # if we're running the analysis on more than one participant...

        # preallocate lists for each task
        groupNull = [[]] * 2
        groupFitData = [[]] * 2

        # for each task...
        for TASK in [0,1]:

            #feedback
            print('starting group level null distribution fits for the ' + taskNames[TASK] + ' task')

            # preallocate sublists for groupNull array
            groupNull[TASK] = [[]] * 5

            # concatenate data from the first two participants to get things started
            groupNull[TASK][0] = np.concatenate((nullCorr[TASK][0],nullCorr[TASK][1]),axis=0)

            # concatenate data from any remaining participants
            if numSubs > 2: # if we're running the analysis on more than two participants...
                for SUB in range(2,numSubs): # for each participant...
                     groupNull[TASK][0] = np.concatenate((groupNull[TASK][0],nullCorr[TASK][SUB]),axis=0) # concatenate

            # fit normal distribution to group null
            groupNull[TASK][1] = [[]] * groupNull[TASK][0].shape[1] # preallocate list for each voxel
            groupNull[TASK][2] = [[]] * groupNull[TASK][0].shape[1] # preallocate list for each voxel
            groupNull[TASK][3] = np.zeros((groupNull[TASK][0].shape[1],), dtype=int) # preallocate vector of zeros (length = voxels) for flagging voxels with bad fits
            for VOX in range(groupNull[TASK][0].shape[1]): # for each voxel... (i.e., each column in groupNull[TASK][0])

                # fit normal distribution - using 3 different methods here as a sanity check
                groupNull[TASK][1][VOX] = [[]] * 2
                groupNull[TASK][1][VOX][0] = np.mean(groupNull[TASK][0][:,VOX])
                groupNull[TASK][1][VOX][1] = np.std(groupNull[TASK][0][:,VOX])
                groupNull[TASK][2][VOX] = stats.kstest(groupNull[TASK][0][:,VOX], "norm", groupNull[TASK][1][VOX]) # get KS goodness of fit

                if groupNull[TASK][2][VOX][1] < alpha:
                    groupNull[TASK][3][VOX] = 1
            groupNull[TASK][4] = pd.DataFrame(data={'numBadFits': [sum(groupNull[TASK][3])], 'propBadFits': [sum(groupNull[TASK][3]) / len(groupNull[TASK][3])]}) # number and proportion of voxels with bad normal fits to the group null distribution

            # make another variable with just the group level fit parameters, goodness of fit results, and summary measures
            groupFitData[TASK] = [[]] * 3
            groupFitData[TASK] = [[]] * 3
            groupFitData[TASK][0] = groupNull[TASK][1] # fit parameters
            groupFitData[TASK][1] = groupNull[TASK][2][VOX] # KS test results
            groupFitData[TASK][2] = groupNull[TASK][4]

            #feedback
            print('*** group level null fits ***')
            print(groupFitData[TASK][2])

    # Get total analysis duration
    endTime = time.time()
    duration = (endTime - startTime) / 60 #[min]
    print('control tasks ISC duration: ' + str(duration))

    # save data
    if saveOutput:
        print('saving output...')
        saveFile = outputFolder + saveName + '.pkl'
        with open(saveFile,'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump([permTest, corrData, groupFitData, duration], f, protocol=4)
        print('output saved here: ' + saveFile)

    print('\n\n\nThe loony bin is closing shop, greetings',
          'from The Combine!\n\n')

#%% GO
if __name__ == '__main__':
    main()