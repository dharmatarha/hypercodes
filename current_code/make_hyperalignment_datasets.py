#!/usr/bin/env conda run -n py27Env python2.7
# -*- coding: utf-8 -*-

"""
Make pymvpa datasets using the pymvpa fmri_dataset() function. Takes EPI nifti files
(generated from generate_post_parseEPI_nifti_files_for_hyperalignment.ipynb for control
tasks and resort_storytelling_time_series.ipynb for storytelling tasks) and whole
brain mask (mni_icbm152_t1_tal_nlin_asym_09c_mask_RESAMPLED.nii) as inputs.
"""

import glob
from mvpa2.suite import *
import mvpa2.base.dataset as mvpads

def main():

    print('\nlets make some datasets\n')

    # define hyperscanning task descriptions
    taskDescrips = ['storytelling_independent',
                    'storytelling_joint',
                    'listening',
                    'reading']

    # input parameters
    debug = False
    task = 2 # see task descriptions above

    if debug:
        dataFolder = '/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/preprocessing/hyperalignment/input_nifti_files/debugging/'
    else:
        dataFolder = '/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/preprocessing/hyperalignment/input_nifti_files/'

    # get list of nifti files in dataFolder
    if task >= 3:
        files = glob.glob(dataFolder + '*storytelling' + str(task) + '*')
    elif task == 1:
        files = glob.glob(dataFolder + '*storytelling*' + '_ind' + '*')
    elif task == 2:
        files = glob.glob(dataFolder + '*storytelling*' + '_joint' + '*')

    # set path to whole brain mask
    maskFile = '/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/misc/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_RESAMPLED.nii'

    # get subject IDs in 'files'
    subIDs = [[]] * len(files)
    for FILE in range(len(files)):
        subIDs[FILE] = files[FILE][files[FILE].find('sub-')+4:files[FILE].find('sub-')+13]

    # initialize full dataset as list
    ds = []

    # load each file with the whole brain mask as a pymvpa dataset and append to list
    for FILE in range(len(files)):
        tmp = fmri_dataset(files[FILE], mask=maskFile)
        tmp.a['subID'] = files[FILE][files[FILE].find('sub-')+4:files[FILE].find('sub-')+13] # append subject ID to dataset
        tmp.a['taskNum'] = int(files[FILE][files[FILE].find('telling')+7]) # append subject ID to dataset
        tmp.a['taskDescription'] = taskDescrips[int(files[FILE][files[FILE].find('telling')+7]) - 1]
        tmp.a['full_dataset_subID_list'] = subIDs
        ds.append(tmp)

    # save as h5 file
    if debug:
        savePath = '/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/preprocessing/hyperalignment/datasets/debug_' + str(taskDescrips[task])
        mvpads.save(ds,savePath)
        print('dataset saved here: ' + savePath)
    else:
        savePath = '/dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/preprocessing/hyperalignment/datasets/' + taskDescrips[task]
        mvpads.save(ds,savePath)
        print('dataset saved here: ' + savePath)

    print('\n\n\nThe loony bin is closing shop, greetings',
          'from The Combine!\n\n')

if __name__ == '__main__':
    main()
