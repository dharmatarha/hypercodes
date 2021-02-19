# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 10:11:35 2018

Interpolates data from each voxel according to the supplied time array.
No extrapolation!

Inputs:
file: path to nifti file
tarray: time array for resampling / interpolation

Output:
saved out nifti file with modified name

@author: adamb
"""


#%% Imports

from nilearn import image as nImage
from nilearn import input_data as nInput
import numpy as np
from scipy import signal
from scipy import interpolate
import os
import argparse


#%% Magic numbers, parameters

def params():
    





return


#%% Load data








def main():

    # parse ipnut arguments
    parser = argparse.ArgumentParser()

    # Input arg "folder"
    parser.add_argument(
        'file',
        type=str,
        help='Nifti file')

    # parse arguments, get list
    args = parser.parse_args()

    # check input
    if not os.path.exists(args.file):
        raise ValueError('Input arg "file" is not found')

    # load img
    img = nImage.load_img(file)
