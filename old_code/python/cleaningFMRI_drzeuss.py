# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 05:53:26 2018

quick and dirty scripting for timeseries extraction and cleaning
using nilearn

@author: adamb
"""

import nilearn
from nilearn import input_data


#%% Declare files we use: EPI + brainmask

# folders for input files, folder for outputs
folder = '/flash/wheatley/adamb/hyperscanning_CBS/' \
         'sub-hid000003_fmriprep/fmriprep/sub-hid000003/ses-pair03/'
outputF = '/flash/wheatley/adamb/hyperscanning_CBS/' \
          'sub-hid000003_test/'

# run 1, MNI space
epi_file = folder + 'func/sub-hid000003_ses-pair03_task-storytelling1_run-01_' \
                    'bold_space-MNI152NLin2009cAsym_preproc.nii.gz'
epi_mask_file = folder + 'func/sub-hid000003_ses-pair03_task-storytelling1_run-01_' \
                         'bold_space-MNI152NLin2009cAsym_brainmask.nii.gz'
# anat, MNI space
anat_file = folder + 'anat/sub-hid000003_ses-pair03_acq-MPRAGE_T1w_space-' \
                     'MNI152NLin2009cAsym_preproc.nii.gz'
anat_mask_file = folder + 'anat/sub-hid000003_ses-pair03_acq-MPRAGE_T1w_space-' \
                          'MNI152NLin2009cAsym_brainmask.nii.gz'

# confounds for run1
confounds_file = outputF + 'run1_confounds_short.csv'

##%% Plot image and mask
#nilearn.plotting.plot_roi(anat_mask_file,
#                          bg_img=anat_file,
#                          cmap='Paired')


#%% Extract timeseries with nilearn NiftiMasker
# Importantly, we don't use signal.clean here as we want to do nuisance
# regression in the same step, and we don't have the confoudns option
# through NiftiMasker

epi_masker = input_data.NiftiMasker(epi_mask_file,
                                            t_r=0.727,
                                            verbose=2)

epi_timeseries = epi_masker.fit_transform(epi_file)


#%% Nuisance regression + standardization

#epi_cleaned = nilearn.signal.clean(epi_file)

















