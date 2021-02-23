#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 10:04:31 2018

USAGE: python3 confoundsTruncate_v2.py [--FWDthreshold] [--motionder] folder 
E.g. python3 confoundsTruncate_v2.py -f 1.5 -m /home/adamb/myFMRIPREPfolder

Script to create a truncated confounds file for nuisance regression

The program lists all standard fmriprep confounds tsv files in the input
folder and uses them as input to a "truncate" function. The "truncate"
part deletes all the columns except:
- FWD
- 6 or 12 motion parameters (depending on --motionder)
- aCompCor first five components
It also creates one extra column for each TR with too much movement (FWD)
More specifically, whenever FWD is above the threshold specified in the
input argument, a separate regressor is included for that specific plus
for the next TR (scrubbing).
These scrubbing columns are named Scrub1, Scrub2,..., ScrubN

Output is the truncated csv file, saved next to the originals

@author: adamb
"""


#%% Imports

import argparse
import os
import glob
import pandas as pd
import numpy as np


#%% Magic numbers, parameters declared here
def params():
    # string used for finding confounds files
    searchTerm = '/**/*confounds.tsv'
    # columns of confounds tables to keep
    columnsToKeep = ['X', 'Y', 'Z', 'RotX', 'RotY', 'RotZ', 'aCompCor00',
                     'aCompCor01', 'aCompCor02', 'aCompCor03',
                     'aCompCor04', 'FramewiseDisplacement', 'CSF']
    # names of motion parameter columns - used later for motion derivatives
    motionCols = ['X', 'Y', 'Z', 'RotX', 'RotY', 'RotZ']
    # add this to names of truncated confoudns tsv files
    saveFileAdd = '_truncated'
    # add this to name for text file listing scrubbed TRs
    saveScrub = '_scrubList.txt'

    return searchTerm, columnsToKeep, motionCols, saveFileAdd, saveScrub


#%% get a list of confound files to work with
def truncate(file, columnsToKeep, motionCols, saveFileAdd,
             FWDthreshold, motionder, saveScrub):

    '''
    Function to create the truncated confounds files

    Inputs:
    file - file path
    columnsToKeep - list of column names that are to be kept
    saveFileAdd - the original file name is appended with this string
    FWDthreshold - threshold value for creating scrubbing columns
    motiondr - flag for 12DOF motion correction factors instead of the
                standard 6 (by adding derivatives)
    saveScrub - file name ending for saving out indices of scrubbed TRs

    Output:
    Truncated confound file is saved as original_filename+saveFileAdd.csv
    '''

    print('\n\n----------------',
          '\nTruncating', file)

    # get folder name, base name, extension name, create savefile name
    fileFolder = os.path.dirname(file)
    fileBase = os.path.splitext(os.path.basename(file))[0]
    savePath = fileFolder + '/' + fileBase + saveFileAdd + '.csv'
    savePathScrub = fileFolder + '/' + fileBase + saveScrub

    # read in file
    table = pd.read_table(file)

    # delete not needed columns
    columns = table.columns.values
    columnsToDrop = [col for col in columns if col not in columnsToKeep]
    table.drop(columnsToDrop, inplace=True, axis=1)

    # Add extra regressors for each TR with FWD above threshold
    #
    # Do this only if we have FramewiseDisplacement among the columns
    # we kept
    if 'FramewiseDisplacement' in columnsToKeep:

        # get indices of all TRs with FWD higher than threshold
        idx = np.where(table.FramewiseDisplacement > FWDthreshold)[0]
        # percentage of TRs affected
        percTR = len(idx)/len(table.FramewiseDisplacement)*100

        # give info to user about number of TRs to be scrubbed
        print('\nFound', str(len(idx)), 'TRs with FWD above the',
              'threshold of', '{0:.2f}'.format(FWDthreshold), 'mm',
              '-', '{0:.2f}'.format(percTR),
              '% of all data')
        # list indices too, save them out into scrub log file
        print('Affected TRs:', str(idx))
        with open(savePathScrub, 'w+') as scrubFile:
            for id in idx:
                print(id, file=scrubFile)

        # go through all indices, create a corresponding regressor,
        # then append the truncated confounds table with it
        for num, i in enumerate(idx):
            # get a basic array of zeros we can alter for all extra regressors
            iRegr = np.zeros(len(table.FramewiseDisplacement))
            # think of the case when the index is that of the last element
            if i != len(table.FramewiseDisplacement):
                iRegr[i:i+2] = np.ones(2)
            else:
                iRegr[i] = 1
            # append the table with the new column
            newColumn = 'Scrub' + str(num)
            table[newColumn] = iRegr

    # if the motion derivatives flag was provided,
    # create the extra confound regressors for them
    # and add them to the confounds file
    if motionder and set(motionCols).issubset(columnsToKeep):
        for colName in motionCols:
            motionConf = np.diff(table[colName], 1)
            newColumn = colName + 'd'
            table[newColumn] = np.append(0, motionConf)
        print('\nAdded difference arrays for the motion parameters to the confounds table')

    # save truncated table
    table.to_csv(savePath, index=False)
    print('\nSaved truncated file to', savePath)

    return


#%% Main
def main():

    # parse ipnut arguments
    parser = argparse.ArgumentParser()

    # Input arg "folder"
    parser.add_argument(
        'folder',
        type=str,
        help='Work folder. The script creates a truncated version of '
             'all fmriprep confounds files found in folder '
             '(search is recursive!)')

    # Input arg "FWDthreshold"
    parser.add_argument(
        '-f',
        '--FWDthreshold',
        type=float,
        nargs='?',
        default=2,
        help='FrameWiseDisplacement threshold for scrubbing in mm.'
             'Default value is 2 mm.')

    # Input arg "motionder" (motion derivatives)
    parser.add_argument(
        '-m',
        '--motionder',
        action='store_true',
        help='Flag for 12DOF motion correction instead of the standard 6')

    # parse arguments, get list
    args = parser.parse_args()

    # check inputs
    if not os.path.exists(args.folder):
        raise ValueError('Input arg "folder" is not found')
    if not (0.1 <= args.FWDthreshold <= 5):
        raise ValueError('Input arg "FWDthreshold" seems very unusual, '
                         'not between the expected range of 0.1 - 5.0 mm')

    # small feedback
    if args.motionder:
        print('\nStarted confoundsTruncate.py with arguments', args.folder,
              'as folder and', args.FWDthreshold, 'as FWDthreshold, ',
              'also including difference arrays for the 6 motion parameters')
    else:
        print('\nStarted confoundsTruncate.py with arguments', args.folder,
              'as folder and', args.FWDthreshold, 'as FWDthreshold')

    # load magic numbers
    searchTerm, columnsToKeep, motionCols, saveFileAdd, saveScrub = params()
    print('\nLoaded preset parameters')

    # get list of counfounds files
    fileList = glob.glob(args.folder + searchTerm, recursive=True)
    print('\nGot list of confound files, found', str(len(fileList)))

    for file in fileList:
        truncate(file,
                 columnsToKeep,
                 motionCols,
                 saveFileAdd,
                 args.FWDthreshold,
                 args.motionder,
                 saveScrub)

    print('\n\n\nThe loony bin is closing shop, greetings',
          'from The Combine!\n\n')

    return


#%% GO
if __name__ == '__main__':
    main()
