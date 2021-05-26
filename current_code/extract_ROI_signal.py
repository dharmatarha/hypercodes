#!/usr/bin/env python3


import numpy as np
from nilearn import image as nImage
from nilearn import datasets as nDatasets
from nilearn.input_data import NiftiLabelsMasker
from matplotlib import pyplot


###########################
# basic params, file paths
###########################

# Base working dir
base_dir = '/home/adamb/PycharmProjects/hypercodes/'

# Atlas file downloaded from templateflow website, already registered to the MNI152NLin2009cAsym template:
# HOCPAL = Harvard-Oxford Coortical Probabilistic Atlas Lateralized
# th25 = thresholding level for deriving a discrete version of the atlas from the probabilistic maps
# dseg = discrete segmentation
atlas_file = base_dir + 'tpl-MNI152NLin2009cAsym_res-02_atlas-HOCPAL_desc-th25_dseg.nii.gz'

# number of expected cortical segments in the atlas defined above:
expected_seg = 48*2  # *2 because it is lateralized

# As there is no label file distributed by templateflow with the atlas, we load the FSL-version of the atlas
# available in nilearn. The name of the atlas for nilearn.datasets.fetch_atlas_harvard_oxford() :
atlas_fsl_id = 'cort-maxprob-thr25-2mm'

# Expected label values for primary auditory cortex, that is, for posterior superior temporal gyrus and Heschl's gyrus, bilaterally
roi_labels = ['Left Superior Temporal Gyrus, posterior division',
              'Right Superior Temporal Gyrus, posterior division',
              'Left Heschl\'s Gyrus (includes H1 and H2)',
              'Right Heschl\'s Gyrus (includes H1 and H2)']
roi_label_values= [19, 20, 89, 90]

# 4D BOLD image for ROI signal extraction:
fmri_image = base_dir + 'sub-hid000002_ses-pair02_task-storytelling3_run-03_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021.nii.gz'


###########################
# load atlases
###########################

# load templateflow atlas img, print some info about it
atlas = nImage.load_img(atlas_file)
print('\n\nLoaded atlas from \n' + atlas_file)
print('\nAffine is: ')
print(atlas.affine)
print('\nImage dimensions:')
print(atlas.shape)

# check the number of regions in the atlas
atlas_data = atlas.get_fdata()
# 0 is the background, don't count that as segmentation label
seg_labels = np.unique(atlas_data)
print('\nThere are ' + str(seg_labels.shape[0]-1) + ' labels (=segments) in the atlas (plus one for background)')
if expected_seg == seg_labels.shape[0]-1:
    print('We seem to be fine, expected ' + str(expected_seg) + '.')
else:
    print('PROBLEM: we expected ' + str(expected_seg) + '!')
    raise ValueError('Wrong number of labels in atlas!')

# load FSL version, get segmentation labels
# IMPORTANT: we want the lateralized version, set "symmetric_split"
fsl_atlas = nDatasets.fetch_atlas_harvard_oxford(atlas_fsl_id, symmetric_split=True, verbose=5)
labels = fsl_atlas['labels']
print('\n\nLoaded FSL atlas ' + atlas_fsl_id)
# check number of labels
print('\nThere are ' + str(len(labels)-1) + ' labels (=segments) in the atlas (plus one for background)')
if expected_seg == len(labels)-1:
    print('We seem to be fine, expected ' + str(expected_seg) + '.')
else:
    print('PROBLEM: we expected ' + str(expected_seg) + '!')
    raise ValueError('Wrong number of labels from FSL atlas!')

# check if the labels of interest are as expected
for i in range(len(roi_labels)):
    # error if ROI labels do not correspond to expected label values in atlas
    if labels.index(roi_labels[i]) != roi_label_values[i]:
        print('Atlas label value (=index) for ' + roi_labels[i] + ' is not the expected value! Investigate!')
print('\nAll ROI labels are as expected, we pair the templateflow atlas with labels from the FSL version.')


###########################
# Extract ROI signals
###########################

# load fmri data
data = nImage.load_img(fmri_image)
print('\nLoaded BOLD image')
print('Affine is: ')
print(data.affine)

# Define a masker object for ROI labels
# The masker would happily do signal cleaning and confound regression,
# but we only use it for signal extraction here as the data has already been preprocessed.
# Important arguments are:
# - "resampling_target" - whether data is to be resampled to atlas or the other way around
# - "strategy" - how the masker should reduce the signals from a region to one time series
masker = NiftiLabelsMasker(labels_img=atlas,
                           background_label=0,
                           standardize=False,
                           resampling_target='data',
                           strategy='mean',
                           verbose=5)

# Extract ROI signals
# Implicitly resamples atlas (labels img) to data dimensions and affine
time_series = masker.fit_transform(data)
time_series_labels = labels[1:]
ts_label_indices = [time_series_labels.index(roi_labels[0]),
                    time_series_labels.index(roi_labels[1]),
                    time_series_labels.index(roi_labels[2]),
                    time_series_labels.index(roi_labels[3]),
                    ]



###########################
# Plotting, look for regions with strange pattern at the beginning
###########################

# Time axis, TRs
trs = [i for i in range(time_series.shape[0])]

# First two ROIs - bilateral STG
pyplot.subplot(1, 2, 1)
pyplot.plot(trs, time_series[:, ts_label_indices[0]])
pyplot.plot(trs, time_series[:, ts_label_indices[1]])
pyplot.title('Preprocessed signal from first two ROIs')
pyplot.xlabel('TRs')
pyplot.ylabel('z-score')
pyplot.legend([roi_labels[0], roi_labels[1]])

# Second two ROIs - bilateral Heschl's
pyplot.subplot(1, 2, 2)
pyplot.plot(trs, time_series[:, ts_label_indices[2]])
pyplot.plot(trs, time_series[:, ts_label_indices[3]])
pyplot.title('Preprocessed signal from last two ROIs')
pyplot.xlabel('TRs')
pyplot.ylabel('z-score')
pyplot.legend([roi_labels[2], roi_labels[3]])
pyplot.show()

