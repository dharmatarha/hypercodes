#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to perform nuisance regression and timeseries extraction on listening
(Run 3) and reading (Run 4) task data.
USAGE: python3 nuisRegr_runs34_2021.py [--mask] folder
The script first lists all functional nii.gz files and truncated confounds files
for the inputted run (3 or 4). Then it applies nilearn's NiftiMasker to each
which calls signal.clean internally, using the supplied confounds file.
Standardization, spatial smoothing and high-pass filtering are performed if
specified so in the magic numbers section.
Current default is high-pass filtering (cutoff at 0.01 Hz)
and spatial smoothing (8 mm fwhm).
Importantly, we only look for run3 and run4 files.
By default, the script uses the brainmask in the fmriprep output
corresponding to the EPI file, change this behavior with the --mask
option, by supplying the mask to be used.
Input arg:
folder: path of data folder
Option --mask, -m:
    A path to the mask to be used OR
    "nilearnMNI" to use the built-in MNI152 brainmask in nilearn
If --mask is not supplied, the script looks for a brainmask in the
standard fmriprep output form and uses that. Be careful with your own
mask, check if its affine matches the EPI file affines
Output is a saved out np array for extracted and cleaned timeseries for
each functional run
Created by adamb on Thu Apr  5 05:53:26 2018
Modified by JDK, March-April 2021
"""

import os
import glob
import argparse
from nilearn import image as nImage
from nilearn import input_data
from nilearn import datasets as nDatasets
import numpy as np
import pandas as pd
import nibabel as nb
from scipy import io as sio


#%% Params, magic numbers
def params():

    # Search terms for glob.glob
    # by default we use the images registered to the MNI 152 nonlinear
    # 2009c template
    epiTerms = ['/**/*run-03_bold_space-'
                'MNI152NLin2009cAsym_preproc.nii.gz',
                '/**/*run-04_bold_space-'
                'MNI152NLin2009cAsym_preproc.nii.gz'
                ]
    brainmaskTerm = 'brainmask.nii.gz'
    confoundsTerm = 'confounds_truncated_2021.csv' #***JD edit*** added the '_2021'
    switchNo = [-1, -2]

    # !! WATCH OUT FOR SETTING THE RIGHT TR VALUE !!
    tr = 0.727
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # glue this to the end of the epi file name for the output file
    saveFileAdd = '_nuisRegr_newMask' #***JD edit*** added the '_2021'
    # switch NaN values (if any) in the confound file to this vale
    NaNValue = 0
    detrending = False
    if detrending:
        saveFileAdd = saveFileAdd + '_detrended'
    # settings for the NiftiMasker object
    masker_verbose = 5 # jd edit: originally set to 2
    masker_cutoff = 0.01  # high-pass filter cutoff [Hz]
    masker_fwhm = 8  # spatial smoothing [mm]
    masker_stand = True  # standardization
    # mni153 nonlin 2009c affine when resampled to 3x3x3
    expectedAffine = np.array([[3, 0, 0, -96],
                               [0, 3, 0, -132],
                               [0, 0, 3, -78],
                               [0, 0, 0, 1]
                               ])
    expectedShape = (65, 77, 65)
    # cutoff when we resample the built-in MNI152 brainmask from nilearn
    # adn have to deal with values between 0 and 1
    maskBoolCutoff = 0.1
    # save the nuisance regressed time series into a mat file as well?
    matSave = True

    return [epiTerms, brainmaskTerm, confoundsTerm, switchNo, tr,
            saveFileAdd, NaNValue, detrending, masker_verbose, masker_cutoff,
            masker_fwhm, masker_stand, expectedAffine, expectedShape,
            maskBoolCutoff, matSave]


#%% Function to list all files in input folder eligible for cleaning
def listFiles(folder, epiTerms, brainmaskTerm, confoundsTerm, switchNo, confoundsFolder):

    '''
    List all eligible files by globbing the given pattern and verifying
    the existence of corresponding confounds file
    '''

    # init final result lists
    epiFiles = []
    maskFiles = []
    confFiles = []
    # go through all given glob prompts
    for term in epiTerms:
        tempList = glob.glob(folder + term, recursive=True)

        # check each result for brainmask + confounds
        for epiFile in tempList:
            fileFolder = os.path.dirname(epiFile)
            fileBase = os.path.splitext(os.path.basename(epiFile))[0]
            # brainmask is only different in the last element of filename
            fileBaseParts = fileBase.split('_')  # split into BIDS parts
            maskFile = '_'.join(fileBaseParts[0:switchNo[0]] +
                                [brainmaskTerm])  # switch last element
            # confounds file is different in last two elements
            confFile = '_'.join(fileBaseParts[0:switchNo[1]] +
                                [confoundsTerm])
            print(confFile)

            # if confound file exists, we append final result lists
            if os.path.exists(confoundsFolder + '/' + confFile):
                epiFiles.append(epiFile)
                confFiles.append(confoundsFolder + '/' + confFile)
                # maskFiles is appended only if file exists
                if os.path.exists(fileFolder + '/' + maskFile):
                    maskFiles.append(fileFolder + '/' + maskFile)
                else:
                    maskFiles.append('')

    return epiFiles, maskFiles, confFiles


#%% function to load MNI brain mask and resample it to affine of our MNI data
def getMask(epiFiles, expectedAffine, expectedShape, maskBoolCutoff):

    '''
    Function to load built-in nilearn MNI template and resample to
    fmriprep MNI output affine (3x3x3 mm and different offsets)
    '''

    # Sanity check:
    # load a sample epi file, get affine, compare with our pre-set one
    epiExample = nImage.load_img(epiFiles[0])
    if not np.array_equal(expectedAffine, epiExample.affine):
        raise ValueError('\nExpected affine does not match affine '
                         'of first EPI file')

    # resample built-in MNI mask to the affine of the fmriprep outputs
    maskMNI = nDatasets.load_mni152_brain_mask()
    maskMNIResampled = nImage.resample_img(maskMNI,
                                           target_affine=expectedAffine,
                                           target_shape=expectedShape)

    # deal with the effects of continuous interpolation, get boolean again
    dataMask = maskMNIResampled.get_data()  # data part of img
    dataMask[dataMask >= maskBoolCutoff] = 1
    dataMask[dataMask < maskBoolCutoff] = 0
    # new img from changed data
    newMaskImg = nb.Nifti1Image(dataMask,
                                maskMNIResampled.affine,
                                maskMNIResampled.header)

    return newMaskImg


#%% Function to deal with NaNs in confounds cvs
def confoundNaN(file, value):

    '''
    Function to swap NaN values in confound csv files to stg else
    '''

    # load csv
    data = pd.read_csv(file)
    # go on if there are NaN values
    if data.isnull().values.any():
        print('\nSwapping NaN values in', file)
        # replace NaN with "value"
        data.fillna(value, inplace=True)
        # save out results with same name
        data.to_csv(file, index=False)

    return


#%% Call nilearn and do masking + cleaning + extraction on epi file
def callNilearn(epiFile, maskImg, confFile, tr, saveFileAdd,
                masker_verbose, detrending, masker_cutoff, masker_fwhm,
                masker_stand, matSave, outFolder):

    '''
    Nuisance regression for given epi input file using a brainmask
    '''

    # create savepath for cleaned timeseries array
    fileFolder = os.path.dirname(epiFile)
    fileBase, ext = os.path.splitext(os.path.basename(epiFile))
    # if file was .nii.gz, we split again to get basename
    if ext == '.gz':
        fileBase = os.path.splitext(os.path.basename(fileBase))[0]
    savePath = outFolder + '/' + fileBase + saveFileAdd + '.nii.gz'

    # Sanity check - compare EPI affine to brainmask affine
    tempImg = nImage.load_img(epiFile)
    if not np.array_equal(tempImg.affine, maskImg.affine):
        print('WARNING! EPI affine does not match resampled MNI mask '
              'affine!\nEPI.affine:', tempImg.affine, '\nMask.affine:',
              maskImg.affine)

    # create a masker object: spatial smoothing, high pass filtering,
    # standardization are all to be performed on masked image
    masker = input_data.NiftiMasker(maskImg,
                                    t_r=tr,
                                    verbose=masker_verbose,
                                    detrend=detrending,
                                    high_pass=masker_cutoff,
                                    smoothing_fwhm=masker_fwhm,
                                    standardize=masker_stand)
    print('\nCreated masker object. NiftiMasker verbosity set to',
          masker_verbose, '\n')

    # clean + extract timeseries
    tseries = masker.fit_transform(epiFile, confounds=confFile)
    print('\nNp array size:', str(tseries.shape))

    # transform back to image, save out
    cleaned_img = masker.inverse_transform(tseries)
    print('\nCleaned img shape:', cleaned_img.header.get_data_shape())
    print('\nCleaned img affine:\n', cleaned_img.affine)

    # # save out img
    # cleaned_img.to_filename(savePath)

    # if mat file output is set to True
    if matSave:
        # save out np array into mat file
        savematFile = outFolder + '/' + fileBase + saveFileAdd + '.mat'
        sio.savemat(savematFile, {'tseries': tseries})

    return savePath


#%% Nuisance regression + standardization

def main():

    # parse ipnut arguments
    parser = argparse.ArgumentParser()

    # Input arg "folder"
    parser.add_argument(
        'folder',
        type=str,
        help='Folder of CBS data. We use glob to list all eligible '
             'epi files before applying nuisance regression and '
             'timeseries extraction. Glob is recursive!')

    # Input arg "outFolder"
    parser.add_argument(
        '-o',
        '--outFolder',
        default='/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/control_tasks/nuisRegr_output_files',
        type=str,
        help='Destination folder for outputs.')

    # Input arg "confoundsFolder"
    parser.add_argument(
        '-c',
        '--confoundsFolder',
        default='/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/control_tasks/nuisRegr_input_files',
        type=str,
        help='Where to look for confounds csv files.')

    # Input option "--mask"
    parser.add_argument(
        '-m',
        '--mask',
        nargs='?',
        type=str,
        default='/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/misc/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask.nii',
        help='Path to brainmask. This mask will be used for the '
             'NiftiMasker object. Default is mni_asym09c_mask_resamp3x3.nii.gz.')

    # parse arguments, get list
    args = parser.parse_args()

    # check inputs
    if not os.path.exists(args.folder):
        raise ValueError('Input arg "folder" is not found')
    if args.mask and not os.path.exists(args.mask):
        raise ValueError('File supplied for --mask is not found')

    # start messages
    print('\n\nnuisRegr.py was started with input folder',
          args.folder)
    if args.mask:
        if args.mask == 'nilearnMNI':
            print('\nWill use nilearn\'s MNI152 mask resampled to 3x3x3')
        else:
            print('\nWill use the brainmask at', args.mask)
    else:
        print('\nWill use fmriprep output brainmask for each EPI file')

    # load magic numbers / preset parameters
    [epiTerms, brainmaskTerm, confoundsTerm,
     switchNo, tr, saveFileAdd, NaNValue,
     detrending, masker_verbose, masker_cutoff,
     masker_fwhm, masker_stand,
     expectedAffine, expectedShape,
     maskBoolCutoff, matSave] = params()
    print('\nLoaded preset parameters. Masker inputs: ',
          '\nmasker_verbose = ' + str(masker_verbose) +
          ';\nmasker_cutoff = ' + str(masker_cutoff) +
          ';\nmasker_detrending = ' + str(detrending) +
          ';\nmasker_fwhm = ' + str(masker_fwhm) +
          ';\nmasker_stand = ' + str(masker_stand))
    if args.mask == 'nilearnMNI':
        print('\nNilearn MNI152 mask affine will be set to\n', expectedAffine)
    if matSave:
        print('\nMasker outputs (nuisance regressed time series) ',
              'will be saved into mat files as well')

    # list eligible files
    epiFiles, maskFiles, confFiles = listFiles(args.folder,
                                               epiTerms,
                                               brainmaskTerm,
                                               confoundsTerm,
                                               switchNo,
                                               args.confoundsFolder)
    print('\nFound', str(len(epiFiles)), 'epi files with '
          'corresponding *' + confoundsTerm + ' file.')

    # get brainmask
    if args.mask == 'nilearnMNI':  # if nilearn MNI mask is used
        # Create a common MNI brain mask with same parameters as our EPI data
        print('\nResampling nilearn\s MNI152 brainmask to match affine '
              'of EPI files. ')
        maskImg = getMask(epiFiles,
                          expectedAffine,
                          expectedShape,
                          maskBoolCutoff)
        print('Resampling done.')
    elif not args.mask:  # if we use fmriprep output brainmasks
        # check if there are missing brainmask files
        idx = [num for num, y in enumerate(maskFiles) if y == '']
        if idx:
            raise ValueError('\n!!!! Problem !!!! At least one EPI file'
                             'has no corresponding brainmask:',
                             [epiFiles[index] for index in idx])
        else:
            print('\nFound brainmask files for all EPI, will use those '
                  'with NiftiMasker')
    else:  # if we point to a brainmask file, open it with nilearn
        try:

            # *** JD 2021 edits ***
            maskImg_OG = nImage.load_img(args.mask)
            maskImg = nImage.resample_img(maskImg_OG,
                                           target_affine=expectedAffine,
                                           target_shape=expectedShape)

            # deal with the effects of continuous interpolation, get boolean again
            dataMask = maskImg.get_data()  # data part of img
            dataMask[dataMask >= maskBoolCutoff] = 1
            dataMask[dataMask < maskBoolCutoff] = 0

            # new img from changed data
            maskImg = nb.Nifti1Image(dataMask,
                                        maskImg.affine,
                                        maskImg.header)

            print('\nLoaded brainmask file', args.mask)
        except:
            raise ValueError('\nCould not load brainmask file:', args.mask)

    # go through confound files, treat NaN values
    print('\nChecking confound files, '
          'replacing NaN values with', str(NaNValue))
    for file in confFiles:
        confoundNaN(file, NaNValue)

    # call the nuisance regr + extraction for each file
    for idx, epiFile in enumerate(epiFiles):
        print('\n\n-------------------------------------')
        print('Calling nilearn\'s NiftiMasker on', epiFile)
        # if we use fmriprep output brainmasks, maskImg will point to
        # the file path of the mask corresponding to current epi, as
        # NiftiMasker can take a nii.gz as input as well
        if not args.mask:
            maskImg = nImage.load_img(maskFiles[idx])
        savePath = callNilearn(epiFile,
                               maskImg,
                               confFiles[idx],
                               tr,
                               saveFileAdd,
                               masker_verbose,
                               detrending,
                               masker_cutoff,
                               masker_fwhm,
                               masker_stand,
                               matSave, args.outFolder)
        print('\nSaved out residual img to', savePath, ',',
              (len(epiFiles)-(idx+1)), 'files left')

    print('\n\n\nThe loony bin is closing shop, greetings',
          'from The Combine!\n\n')


#%% GO
if __name__ == '__main__':
    main()
