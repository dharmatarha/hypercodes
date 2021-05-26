import os
import numpy as np
from nilearn import image as nImage
from nilearn.input_data import NiftiMasker
from nilearn import plotting
from matplotlib import pyplot


'''
Script to check for "bad" voxels in a given functional, 4D nifti image. 
Just edit the file locations and the subject id in the first section and then you should be good to go.

The script expects a preprocessed nifti image (after nuisance regression, standardization, etc.) registered to the 
MNI152NLin2009cAsym template (default fmriprep). Thus, for time series extraction and for plotting it requires paths 
to the brain mask and to the anatomical image for that template, both available from 
https://www.templateflow.org/browse/

The script checks for bad voxels in three simple, but different ways:
- STD method: A voxel is "bad" if its std in the first "N" TRs is larger than 
its (std in the second half of the run) * (some threshold)
- Min-max range method: A voxel is "bad" if its minimum-maximum range in the first "N" TRs is larger than 
its (min-max range in the second half of the run) * (some threshold)
- Mean method: A voxel is "bad" if the absolute difference between its mean in the first "N" TRs and 
its mean in the second half of the run is larger than some threshold

The parameters for voxel detection are set in the "Detect "bad" voxels" section. 
N = 30 and thresholds = (2, 1.4, 0.6) seemed to provide reasonable detections.

The script prints the number of detected "bad" voxels according to each method into standard out. 
Also prints the correlations between the different methods (correlations of binary vectors for all 
voxels encoding "bad" and "good").

The script saves out the following files into "output_dir" (one output folder is created for each subject, 
based on "sub_id"):
- Example time series. By default, it plots 60 (=plot_no) time series in groups of 6 (=subplot_no), 
saved into .png files, separately for each detection method. Useful for tuning the detection parameters 
and checking the kind of "badness" we can detect with a given method.
- Brain masks of "bad" voxels, according to each method. Useful for loading it onto the anatomical template in any 
brain viewer and seeing the extent and location of detected voxels.
- Sample brain images for "bad" voxels. Using nilearn.plotting, we generate three standard slices on the 
anatomical image, centered on the largest cluster of detected voxels (default nilearn behavior).
'''


###############################
# File locations
###############################

# Base working dir
base_dir = '/home/adamb/PycharmProjects/hypercodes/'
# brain mask from templateflow
brain_mask_file = base_dir + 'tpl-MNI152NLin2009cAsym_res-01_desc-brain_mask.nii.gz'
# anatomical template for background, from templateflow
anat_file = base_dir + 'tpl-MNI152NLin2009cAsym_res-01_T1w.nii.gz'
# subject ID
sub_id = 'hid003'
# fmri image
fmri_img_file = base_dir + 'sub-hid000003_ses-pair03_task-storytelling3_run-03_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021.nii.gz'
# output folder for subject
output_dir = base_dir + sub_id + '/'

# make output folder if it does not exist
if not os.path.exists(output_dir):
    os.mkdir(output_dir)



###############################
# Load data, mask it,
# extract time series
###############################

# load image and mask
brain_mask = nImage.load_img(brain_mask_file)
data = nImage.load_img(fmri_img_file)

# extract time series, don't do any more signal cleaning
masker = NiftiMasker(brain_mask,
                     standardize=False,
                     verbose=5,
                     target_affine=data.affine,
                     target_shape=data.shape[0:3])
time_series = masker.fit_transform(data)



###############################
# Detect "bad" voxels
###############################
# Methods:
# - check each voxel for large variance in first "N" TRs, relative to second half of TRs
# - check each voxel for a large min-max range in first "N" TRs, relative to second half of TRs
# - check each voxel for a different mean in first "N" TRs, relative to second half of TRs

# params
# number of TRs in time series
tr_no = time_series.shape[0]
# number of TRs in the supposedly "bad" part of the time series in the beginning
N = 30
# detection threshold for STD method: voxels where STD is larger in the first "N" TRs than
# threshold*STD in the second half, are considered "bad"
std_threshold = 2
# detection threshold for Min-Max method: if min-man range is threshold times larger in
# first "N" TRs than in the second half, voxel is "bad"
minmax_threshold = 1.4
# detection threshold for Mean-method: if the difference between the absolute mean in first "N" TRs
# and the absolute mean in the second half is larger than threshold, voxel is "bad"
mean_threshold = 0.8

# STD method:
first_segment_std = np.std(time_series[:N, :], axis=0)
second_segment_std = np.std(time_series[round(tr_no/2):, :], axis=0)
bad_voxels_std = first_segment_std > (second_segment_std*std_threshold)
print('\nFound ' + str(np.sum(bad_voxels_std)) +
      ' "bad" voxels with STD method (N = '
      + str(N) + ' and threshold = ' + str(std_threshold) + ')')
bad_voxels_std = bad_voxels_std.astype(float)  # into float from boolean

# Min-max method:
first_segment_minmax = time_series[:N, :].max(axis=0) - time_series[:N, :].min(axis=0)
second_segment_minmax = time_series[round(tr_no/2):, :].max(axis=0) - time_series[round(tr_no/2):, :].min(axis=0)
bad_voxels_minmax = first_segment_minmax > (second_segment_minmax*minmax_threshold)
print('\nFound ' + str(np.sum(bad_voxels_minmax)) +
      ' "bad" voxels with Min-Max range method (N = '
      + str(N) + ' and threshold = ' + str(minmax_threshold) + ')')
bad_voxels_minmax = bad_voxels_minmax.astype(float)  # into float from boolean

# Mean method:
first_segment_mean = np.mean(time_series[:N, :], axis=0)
second_segment_mean = np.mean(time_series[round(tr_no/2):, :], axis=0)
bad_voxels_mean = np.abs(first_segment_mean-second_segment_mean) > mean_threshold
print('\nFound ' + str(np.sum(bad_voxels_mean)) +
      ' "bad" voxels with Mean difference method (N = '
      + str(N) + ' and threshold = ' + str(mean_threshold) + ')')
bad_voxels_mean = bad_voxels_mean.astype(float)  # into float from boolean

# Get correlations across different methods
tmp = np.stack((bad_voxels_std, bad_voxels_minmax, bad_voxels_mean), axis=1)
corr_matrix = np.corrcoef(tmp.T)
print('\n\nCorrelations across detection methods:')
print(corr_matrix)



###############################
# Plot example time series
###############################
# Save out "plot_no" examples into .png files in "base_dir"
# Do it separately for each detection method

# number of example time series to plot for each detection method
# MUST BE MULTIPLE OF "SUBPLOT_NO!
plot_no = 60
# number of subplots in one plot
subplot_no = 6
# Time axis
TRs = [i for i in range(tr_no)]

# loop through detection methods
for detection_method in ['std', 'minmax', 'mean']:
    if detection_method == 'std':
        voxel_indices = np.where(bad_voxels_std == 1)[0]  # indices of "bad" voxels
    elif detection_method == 'minmax':
        voxel_indices = np.where(bad_voxels_minmax == 1)[0]  # indices of "bad" voxels
    elif detection_method == 'mean':
        voxel_indices = np.where(bad_voxels_mean == 1)[0]  # indices of "bad" voxels

    # get "plot_no" indices randomly from "voxel_indices"
    indices_to_plot = voxel_indices[np.random.choice(voxel_indices.shape[0], plot_no, replace=False)]

    # plot time series with subplots
    for plot_file in range(int(plot_no/subplot_no)):
        tmp_indices = indices_to_plot[plot_file*subplot_no:(plot_file+1)*subplot_no]
        save_file_png = output_dir + 'sub_' + sub_id + '_samples_' + detection_method + '_' + str(plot_file) + '.png'
        for i in range(subplot_no):
            pyplot.subplot(int(subplot_no/2), 2, i+1)
            pyplot.plot(TRs, time_series[:, tmp_indices[i]])
            pyplot.xlabel('TRs')
            pyplot.ylabel('z-score')
        pyplot.savefig(save_file_png)
        pyplot.close()



###############################
# Get masks showing the "bad" voxels
# Save Nifti images
###############################

# STD method:
# get a brain image from "bad_voxels", save into a ".nii"
bad_voxel_mask_std = masker.inverse_transform(bad_voxels_std)
# save path for nifti containing a mask of bad voxels
save_file = output_dir + 'sub_' + sub_id + '_bad_voxels_std_n' + str(N) + '_th' + str(std_threshold) + '.nii'
bad_voxel_mask_std.to_filename(save_file)

# Min-max range method:
# get a brain image from "bad_voxels", save into a ".nii"
bad_voxel_mask_minmax = masker.inverse_transform(bad_voxels_minmax)
# save path for nifti containing a mask of bad voxels
save_file = output_dir + 'sub_' + sub_id + '_bad_voxels_minmax_n' + str(N) + '_th' + str(minmax_threshold) + '.nii'
bad_voxel_mask_minmax.to_filename(save_file)

# Min-max range method:
# get a brain image from "bad_voxels", save into a ".nii"
bad_voxel_mask_mean = masker.inverse_transform(bad_voxels_mean)
# save path for nifti containing a mask of bad voxels
save_file = output_dir + 'sub_' + sub_id + '_bad_voxels_mean_n' + str(N) + '_th' + str(mean_threshold) + '.nii'
bad_voxel_mask_mean.to_filename(save_file)



###############################
# Plot masks
###############################

# plot mask from STD method, save into file
standard_fig = pyplot.figure(figsize=(10, 5))  # get figure first, so we can set size
figure_std = plotting.plot_roi(bad_voxel_mask_std, bg_img=anat_file, figure=standard_fig)
figure_std.title('"Bad" voxels detected by STD method, subject ' + sub_id)
figure_std.savefig(output_dir + 'brain_slices_STD_method_sub_' + sub_id + '.png')
figure_std.close()

# plot mask from Min-Max method, save into file
standard_fig = pyplot.figure(figsize=(10, 5))  # get figure first, so we can set size
figure_minmax = plotting.plot_roi(bad_voxel_mask_minmax, bg_img=anat_file, figure=standard_fig)
figure_minmax.title('"Bad" voxels detected by Min-Max range method, subject ' + sub_id)
figure_minmax.savefig(output_dir + 'brain_slices_minmax_method_sub_' + sub_id + '.png')
figure_minmax.close()

# plot mask from Mean method, save into file
standard_fig = pyplot.figure(figsize=(10, 5))  # get figure first, so we can set size
figure_mean = plotting.plot_roi(bad_voxel_mask_mean, bg_img=anat_file, figure=standard_fig)
figure_mean.title('"Bad" voxels detected by Mean method, subject ' + sub_id)
figure_mean.savefig(output_dir + 'brain_slices_mean_method_sub_' + sub_id + '.png')
figure_mean.close()
