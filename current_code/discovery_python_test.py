#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is to test out running a python script on Discovery.
The goal is to load data from the scratch directory. Do something brief
to confirm that we've loaded that data. Then save some confirmation file
in JD's home directory on Discovery.
"""

import scipy.io as sio
import pickle

def main():

    # set data folder
    folder = '/dartfs-hpc/rc/home/z/f00589z/'
    # folder = '/dartfs-hpc/scratch/f00589z/control_timeseries/'
    file = 'sub-hid000002_ses-pair02_task-storytelling3_run-03_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021_interp.mat'
    outFolder = '/dartfs-hpc/rc/home/z/f00589z/control_ISC_output'


    # load data
    tmp = sio.loadmat(folder + file) # load .mat file
    data = tmp['tseries'] # get timeseries data
    printString = 'loaded ' + str(data.shape[0]) + ' x ' + str(data.shape[1]) + ' timeseries'
    print(printString)

    # save something just to make sure we can
    saveFile = folder + 'control_ISC_output/discovery_python_test_saveFile.pkl'
    with open(saveFile, 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([printString], f)

#%% GO
if __name__ == '__main__':
    main()