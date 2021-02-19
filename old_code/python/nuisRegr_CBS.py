# -*- coding: utf-8 -*-
"""
Script to perform nuisance regression and timeseries extraction on CBS data

USAGE: python3 nuisRegr_CBS.py folder

The script first lists all functional runs for which we have EPI data,
extracted brainmask and truncated confounds files. Then it applies
nilearn's NiftiMasker to each which calls signal.clean internally, using
the supplied confounds file.
Importantly, we only look for run1 and run2 files.

Input:
folder - path of CBS data folder

Output is a saved out np array for extracted and cleaned timeseries for
each functional run


Created on Thu Apr  5 05:53:26 2018
@author: adamb
"""

import argparse
import nilearn
from nilearn import input_data
import os
import glob
import numpy as np
import pandas as pd


#%% Params, magic numbers
def params():
    searchTerms = ['/**/*run-01_bold_space-'
                   'MNI152NLin2009cAsym_preproc.nii.gz',
                   '/**/*run-02_bold_space-'
                   'MNI152NLin2009cAsym_preproc.nii.gz'
                   ]
    brainMaskTerm = 'brainmask.nii.gz'
    confoundsTerm = 'confounds_truncated.csv'
    switchNo = [-1, -2]
    tr = 0.727
    saveFileAdd = '_tseries'
    NaNValue = 0
    masker_verbose = 2
    masker_cutoff = 0.01
    masker_fwhm = 8
    masker_stand = True

    return [searchTerms, brainMaskTerm, confoundsTerm, switchNo, tr,
            saveFileAdd, NaNValue, masker_verbose, masker_cutoff,
            masker_fwhm, masker_stand]


#%% Function to list all files in input folder eligible for cleaning
def listFiles(folder, searchTerms, brainMaskTerm, confoundsTerm, switchNo):

    '''
    List all eligible files by globbing the given pattern and verifying
    the existence of corresponding brainmask and confounds file
    '''

    # init final result lists
    epiFiles = []
    maskFiles = []
    confFiles = []
    # go through all given glob prompts
    for term in searchTerms:
        tempList = glob.glob(folder + term, recursive=True)

        # check each result for brainmask + confounds
        for epiFile in tempList:
            fileFolder = os.path.dirname(epiFile)
            fileBase = os.path.splitext(os.path.basename(epiFile))[0]
            # brainmask is only different in the last element of filename
            fileBaseParts = fileBase.split('_')  # split into BIDS parts
            maskFile = '_'.join(fileBaseParts[0:switchNo[0]] +
                                [brainMaskTerm])  # switch last element
            # confounds file is different in last two elements
            confFile = '_'.join(fileBaseParts[0:switchNo[1]] +
                                [confoundsTerm])

            # if both exist, we append final result lists
            if os.path.exists(fileFolder + '/' + maskFile):
                if os.path.exists(fileFolder + '/' + confFile):
                    epiFiles.append(epiFile)
                    maskFiles.append(fileFolder + '/' + maskFile)
                    confFiles.append(fileFolder + '/' + confFile)

    return epiFiles, maskFiles, confFiles


#%% Function to deal with NaNs in confounds cvs
def confoundNaN(file, value):

    '''
    Function to swap NaN values in confound csv files to stg else
    '''

    # load csv
    data = pd.read_csv(file)
    # go on if there are NaN values
    if data.isnull().values.any():
        print('\nFound NaN in', file)
        # replace NaN with "value"
        data.fillna(value)
        # save out results with same name
        data.to_csv(file, index=False)

    return


#%% Call nilearn and do masking + cleaning + extraction on epi file
def callNilearn(epiFile, maskFile, confFile, tr, saveFileAdd,
                masker_verbose, masker_cutoff, masker_fwhm,
                masker_stand):

    '''
    Nuisance regression + timeseries extraction for given epi input file
    '''

    # create savepath for cleaned timseries array
    fileFolder = os.path.dirname(epiFile)
    fileBase = os.path.splitext(os.path.basename(epiFile))[0]
    savePath = fileFolder + '/' + fileBase + saveFileAdd + '.nii.gz'

    # create the masker object: spatial smoothing, high pass filtering,
    # standardization
    masker = input_data.NiftiMasker(maskFile,
                                    t_r=tr,
                                    verbose=masker_verbose,
                                    high_pass=masker_cutoff,
                                    smoothing_fwhm=masker_fwhm,
                                    standardize=masker_stand)
    print('Created masker object')

    # clean + extract timeseries
    tseries = masker.fit_transform(epiFile, confounds=confFile)
    print('\nNp array size:', str(tseries.shape))

    # transform back to image, save out
    cleaned_img = masker.inverse_transform(tseries)
    print('Cleaned img shape:', cleaned_img.get_data_shape())
    print('Img affine:', cleaned_img.affine)

    # save out img
    cleaned_img.to_filename(savePath)

    # save out np array
    np.save(savePath, tseries)

    return savePath


#%% Nuisance regression + standardization

def main():

    # parse ipnut arguments
    parser = argparse.ArgumentParser()

    # Input arg "folder"
    parser.add_argument(
        'folder',
        type=str,
        help='Folder of CBS data. We use glob to list all eligible'
             'epi files before applying nuisance regression and'
             'timeseries extraction. Glob is recursive!')

    # parse arguments, get list
    args = parser.parse_args()

    # check inputs
    if not os.path.exists(args.folder):
        raise ValueError('Input arg "folder" is not found')

    # start message
    print('\nnuisRegr_CBS was started with input folder =',
          args.folder)

    # load magic numbers / preset parameters
    [searchTerms, brainMaskTerm, confoundsTerm,
     switchNo, tr, saveFileAdd, NaNValue,
     masker_verbose, masker_cutoff,
     masker_fwhm, masker_stand] = params()
    print('\nLoaded parameters, magic numbers. Masker inputs: ',
          'masker_verbose = ' + str(masker_verbose) +
          '; masker_cutoff = ' + str(masker_cutoff) +
          '; masker_fwhm = ' + str(masker_fwhm) +
          '; masker_stand = ' + str(masker_stand))

    # list eligible files
    epiFiles, maskFiles, confFiles = listFiles(args.folder,
                                               searchTerms,
                                               brainMaskTerm,
                                               confoundsTerm,
                                               switchNo)
    print('\nFound', str(len(epiFiles)), 'epi files with corresponding '
          'mask and confounds file. Starting nuisance regression + '
          'timeseries extraction on them, one-by-one...')

    # go through confound files, treat NaN values
    print('\nChecking confound files for NaN values, '
          'replacing them with', str(NaNValue))
    for file in confFiles:
        confoundNaN(file, NaNValue)

    # call the nuisance regr + extraction for each file
    for idx, epiFile in enumerate(epiFiles):
        print('\n-------------------------------------')
        print('\nCalling milearn\'s NiftiMasker/signal.clean on', epiFile)
        savePath = callNilearn(epiFile,
                               maskFiles[idx],
                               confFiles[idx],
                               tr,
                               saveFileAdd,
                               masker_verbose,
                               masker_cutoff,
                               masker_fwhm,
                               masker_stand)
        print('Saved out np array to', savePath, ',',
              (len(epiFiles)-(idx+1)), 'files left')

    print('\n\n\nThe loony bin is closing shop, greetings',
          'from The Combine!\n\n')


#%% GO
if __name__ == '__main__':
    main()
