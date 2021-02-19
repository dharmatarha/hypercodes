# -*- coding: utf-8 -*-
"""
Created on Wed May  2 06:37:40 2018

@author: adamb
"""


#%% Imports

import os
from nilearn import input_data as nInput
from nilearn import plotting as nPlotting
from nilearn import image as nImage
from nilearn import datasets
from nilearn import surface
from nilearn import plotting
from scipy import io as sio
import numpy as np


fsaverage = datasets.fetch_surf_fsaverage5()


#%% 
def surfPlot(file):

    # load img
    folder = '/media/adamb/TOSHIBA_EXT/testConfounds/results_First/'
    img = nImage.load_img(folder + file)

    output_img = '/home/adamb/Desktop/surfOutput.png'

    texture = surface.vol_to_surf(img,
                                  fsaverage.pial_left,
                                  interpolation='linear')

    plotting.plot_surf_stat_map(fsaverage.infl_left,
                                texture,
                                hemi='left',
                                threshold=0.01,
                                bg_map=fsaverage.sulc_right,
                                output_file=output_img)

    return


#%%
def plotOneImg(file):

    #%% File handling

    # get to right pwd
    folder = '/home/adamb/Desktop/testConfounds/'
    os.chdir(folder)
    # specify background img, brainmask img, tr, verbosity
    backg_img = '/home/adamb/Desktop/testConfounds/mni_icbm152_t1_tal_nlin_asym_09c.nii'
    brainmask_img = '/home/adamb/Desktop/testConfounds/mni_asym09c_mask_resamp3x3.nii.gz'

    imgFile = folder + file

    # plotting
    savePlotF = folder + 'currentImg_2.png'
    nPlotting.plot_roi(roi_img=imgFile,
                            bg_img=backg_img,
#                            threshold=1,
                            colorbar=False,
                            display_mode='ortho',
                            cut_coords=[-3, 10, 56],  # [-62, -59, -47, -2]
                            annotate=True,
#                            title=plotTitle,
                            output_file=savePlotF)

    return



def plotNiiSlices(pairN, runN):

    #%% File handling

    # get to right pwd
    folder = '/home/adamb/Desktop/testConfounds/results_Second/'
    os.chdir(folder)
    
    # specify background img, brainmask img, tr, verbosity
    backg_img = '/home/adamb/Desktop/testConfounds/mni_icbm152_t1_tal_nlin_asym_09c.nii'
    brainmask_img = '/home/adamb/Desktop/testConfounds/mni_asym09c_mask_resamp3x3.nii.gz'

    # exact mat files
    files = [folder + '/pair' + str(pairN) + '/RsqFDR_run' + str(runN) +
             '_1_img.nii.gz',
             folder + '/pair' + str(pairN) + '/RsqFDR_run' + str(runN) +
             '_2_img.nii.gz']

    # plotting
    savePlotF = folder + '/pair' + str(pairN) + '/RsqFDR_run' + str(runN) +\
                 '_1_img.png'
    nPlotting.plot_stat_map(stat_map_img=files[0],
                            bg_img=backg_img,
#                            threshold=0.04,
                            colorbar=True,
                            display_mode='x',
                            cut_coords=[-62, -2, 3, 47],
                            annotate=True,
#                            title=plotTitle,
                            output_file=savePlotF)


    savePlotF = folder + '/pair' + str(pairN) + '/RsqFDR_run' + str(runN) +\
                 '_2_img.png'
    nPlotting.plot_stat_map(stat_map_img=files[1],
                            bg_img=backg_img,
#                            threshold=0.04,
                            colorbar=True,
                            display_mode='x',
                            cut_coords=[-62, -2, 3, 47],
                            annotate=True,
#                            title=plotTitle,
                            output_file=savePlotF)


    return


#%% All in all
def plotRsq(pairN, runN):

    #%% File handling

    # get to right pwd
    folder = '/home/adamb/Desktop/testConfounds'
    os.chdir(folder)
    
    # specify background img, brainmask img, tr, verbosity
    backg_img = folder + '/mni_icbm152_t1_tal_nlin_asym_09c.nii'
    brainmask_img = folder + '/mni_asym09c_mask_resamp3x3.nii.gz'
    fit_img = folder + '/results_First/singleEPIMNIsample.nii'
    tr = 1.9
    masker_verbose = 2

    ## specify our target story
    #pairN = 2
    #runN = 1

    # exact mat files
    files = [folder + '/results_Second/pair' + str(pairN) + '/coupling_run' + str(runN) +
             '_1_FDR.mat',
             folder + '/results_Second/pair' + str(pairN) + '/coupling_run' + str(runN) +
             '_2_FDR.mat']

    # load data
    part1 = sio.loadmat(files[0])
    part2 = sio.loadmat(files[1])


    #%% Prepare masker

    # create masker object using the 3x3x3 mm MNINlin mask
    masker = nInput.NiftiMasker(brainmask_img,
                                t_r=tr,
                                verbose=masker_verbose)
    # we need to fit the anatomical to it so we can use it with other data as well
    sampleImg = masker.fit(fit_img)


    #%% Fit what we want and plot it

    # part 1
    RsqFDR_1 = part1['RsqFDR']
    RsqFDR_1_img = masker.inverse_transform(RsqFDR_1)
    saveF = folder + '/results_Second/pair' + str(pairN) + '/RsqFDR_run' + str(runN) +\
        '_1_img.nii.gz'
    RsqFDR_1_img.to_filename(saveF)
    plotTitle = 'Pair ' + str(pairN) + ', run ' + str(runN) + ', part 1'
    savePlotF = folder + '/results_Second/pair' + str(pairN) + '/RsqFDR_run' + str(runN) +\
                 '_1_img.png'
    nPlotting.plot_stat_map(stat_map_img=saveF,
                            bg_img=backg_img,
#                            threshold=0.04,
                            colorbar=True,
                            display_mode='ortho',
                            cut_coords=[55, -5, 10],
                            annotate=True,
                            title=plotTitle,
                            output_file=savePlotF)

    # part 2
    RsqFDR_2 = part2['RsqFDR']
    RsqFDR_2_img = masker.inverse_transform(RsqFDR_2)
    saveF = folder + '/results_Second/pair' + str(pairN) + '/RsqFDR_run' + str(runN) +\
             '_2_img.nii.gz'
    RsqFDR_2_img.to_filename(saveF)
    plotTitle = 'Pair ' + str(pairN) + ', run ' + str(runN) + ', part 2'
    savePlotF = folder + '/results_Second/pair' + str(pairN) + '/RsqFDR_run' + str(runN) +\
                 '_2_img.png'
    nPlotting.plot_stat_map(stat_map_img=saveF,
                            bg_img=backg_img,
#                            threshold=0.04,
                            colorbar=True,
                            display_mode='ortho',
                            cut_coords=[55, -5, 10],
                            annotate=True,
                            title=plotTitle,
                            output_file=savePlotF)

    return


#%% Group images

def plotMeanRsq():

    #%% File handling

    # get to right pwd
    folder = '/home/adamb/Desktop/testConfounds'
    os.chdir(folder)

    # specify background img, brainmask img, tr, verbosity
    backg_img = folder + '/mni_icbm152_t1_tal_nlin_asym_09c.nii'
    brainmask_img = folder + '/mni_asym09c_mask_resamp3x3.nii.gz'
    tr = 1.9
    masker_verbose = 2

    ## specify our target story
    #pairN = 2
    #runN = 1

    # exact mat files
    file = folder + '/MeanMaps.mat'

    # load data
    part = sio.loadmat(file)


    #%% Prepare masker

    # create masker object using the 1x1x1 mm MNINlin mask
    masker = nInput.NiftiMasker(brainmask_img,
                                t_r=tr,
                                verbose=masker_verbose)
    # we need to fit the anatomical to it so we can use it with other data as well
    anatImg = masker.fit(backg_img)


    #%% Fit what we want and plot it

    # Joint stories, group
    RsqMeanJoint = part['RsqMeanJoint']
    RsqMeanJoint_img = masker.inverse_transform(RsqMeanJoint)
    saveF = folder + '/RsqMeanJoint_img.nii'
    RsqMeanJoint_img.to_filename(saveF)
    plotTitle = 'Rsquared, group, joint stories'
    savePlotF = folder + '/RsqMeanJoint_img.png'
    nPlotting.plot_stat_map(stat_map_img=saveF,
                            bg_img=backg_img,
                            threshold=0.07,
                            colorbar=True,
                            display_mode='x',
                            cut_coords=6,
                            annotate=True,
                            title=plotTitle,
                            output_file=savePlotF)

    # Joint stories, group
    RsqMeanInd = part['RsqMeanInd']
    RsqMeanInd_img = masker.inverse_transform(RsqMeanInd)
    saveF = folder + '/RsqMeanInd_img.nii'
    RsqMeanInd_img.to_filename(saveF)
    plotTitle = 'Rsquared, group, ind stories'
    savePlotF = folder + '/RsqMeanInd_img.png'
    nPlotting.plot_stat_map(stat_map_img=saveF,
                            bg_img=backg_img,
                            threshold=0.07,
                            colorbar=True,
                            display_mode='x',
                            cut_coords=6,
                            annotate=True,
                            title=plotTitle,
                            output_file=savePlotF)

    return


def plotMeanRsqNonFDR():

    #%% File handling

    # get to right pwd
    folder = '/home/adamb/Desktop/testConfounds'
    os.chdir(folder)

    # specify background img, brainmask img, tr, verbosity
    backg_img = folder + '/mni_icbm152_t1_tal_nlin_asym_09c.nii'
    brainmask_img = folder + '/mni_asym09c_mask_resamp3x3.nii.gz'
    fit_img = folder + '/results_First/singleEPIMNIsample.nii'
    tr = 1.9
    masker_verbose = 2

    ## specify our target story
    #pairN = 2
    #runN = 1

    # exact mat files
#    file = folder + '/results_First/FDRImages.mat'
    file = folder + '/voxelLagsFDR.mat'

    # load data
    part = sio.loadmat(file)


    #%% Prepare masker

    # create masker object using the 3x3x3 mm MNINlin mask
    masker = nInput.NiftiMasker(brainmask_img,
                                t_r=tr,
                                verbose=masker_verbose)
    # we need to fit the anatomical to it so we can use it with other data as well
    anatImg = masker.fit(fit_img)


    #%% Fit what we want and plot it

    # All, group
    voxelEarlyFDR = part['mask']
#    RsqMean = np.power(RsqMean, 0.5)
    voxelEarlyFDR_img = masker.inverse_transform(voxelEarlyFDR)
    saveF = folder + '/voxelLagMask_img.nii'
    voxelEarlyFDR_img.to_filename(saveF)

#    # Joint, group
#    RsqMeanJointFDR = part['RsqMeanJointFDR']
##    RsqMean = np.power(RsqMean, 0.5)
#    RsqMeanJointFDR_img = masker.inverse_transform(RsqMeanJointFDR)
#    saveF = folder + '/results_First/RsqMeanJointFDR_img.nii'
#    RsqMeanJointFDR_img.to_filename(saveF)
##
#    # Ind, group
#    RsqMeanIndFDR = part['RsqMeanIndFDR']
##    RsqMean = np.power(RsqMean, 0.5)
#    RsqMeanIndFDR_img = masker.inverse_transform(RsqMeanIndFDR)
#    saveF = folder + '/results_First/RsqMeanIndFDR_img.nii'
#    RsqMeanIndFDR_img.to_filename(saveF)

    return

#    plotTitle = 'Corr coeffs, group, joint stories, right'
#    savePlotF = folder + '/RsqMeanJoint_img_real_right.svg'
#    nPlotting.plot_stat_map(stat_map_img=saveF,
#                            bg_img=backg_img,
##                            threshold=0.25,
#                            colorbar=True,
#                            display_mode='ortho',
#                            cut_coords=[55,-5,10],
#                            annotate=True,
#                            title=plotTitle,
#                            output_file=savePlotF)

#    # Ind stories, group
#    RsqMeanInd = part['RsqMeanInd']
#    RsqMeanInd = np.power(RsqMeanInd, 0.5)
#    RsqMeanInd_img = masker.inverse_transform(RsqMeanInd)
#    saveF = folder + '/RsqMeanInd_img_real.nii'
#    RsqMeanInd_img.to_filename(saveF)
#    plotTitle = 'Corr coeffs, group, ind stories, right'
#    savePlotF = folder + '/RsqMeanInd_img_real_right.svg'
#    nPlotting.plot_stat_map(stat_map_img=saveF,
#                            bg_img=backg_img,
#                            threshold=0.25,
#                            colorbar=True,
#                            display_mode='ortho',
#                            cut_coords=[55,-5,10],
#                            annotate=True,
#                            title=plotTitle,
#                            output_file=savePlotF)
#
#    # Joint - Ind stories, Contrast
#    RsqMeanContrast = RsqMeanJoint-RsqMeanInd
#    RsqMeanContrast_img = masker.inverse_transform(RsqMeanContrast)
#    saveF = folder + '/RsqMeanContrast_img_real.nii'
#    RsqMeanContrast_img.to_filename(saveF)
#    plotTitle = 'Corr coeffs, joint-ind contrast, medial'
#    savePlotF = folder + '/RsqMeanContrast_img_real_medial.svg'
#    nPlotting.plot_stat_map(stat_map_img=saveF,
#                            bg_img=backg_img,
#                            threshold=0.20,
#                            colorbar=True,
#                            display_mode='ortho',
#                            cut_coords=[0,-5,10],
#                            annotate=True,
#                            title=plotTitle,
#                            output_file=savePlotF)

#    return


def plotMeanRsqContrast():

    #%% File handling

    # get to right pwd
    folder = '/home/adamb/Desktop/testConfounds'
    os.chdir(folder)

    # specify background img, brainmask img, tr, verbosity
    backg_img = folder + '/mni_icbm152_t1_tal_nlin_asym_09c.nii'
    brainmask_img = folder + '/mni_asym09c_mask_resamp3x3.nii.gz'
    tr = 1.9
    masker_verbose = 2

    ## specify our target story
    #pairN = 2
    #runN = 1

    # exact mat files
    file = folder + '/MeanMapsNonFDR.mat'

    # load data
    part = sio.loadmat(file)


    #%% Prepare masker

    # create masker object using the 1x1x1 mm MNINlin mask
    masker = nInput.NiftiMasker(brainmask_img,
                                t_r=tr,
                                verbose=masker_verbose)
    # we need to fit the anatomical to it so we can use it with other data as well
    anatImg = masker.fit(backg_img)


    #%% Fit what we want and plot it

    # Joint stories, group
    RsqMeanContrast1 = part['RsqMeanContrast1']
    RsqMeanContrast1 = np.power(RsqMeanContrast1, 0.5)
    RsqMeanContrast1_img = masker.inverse_transform(RsqMeanContrast1)
    saveF = folder + '/RsqMeanContrast1_img_real.nii'
    RsqMeanContrast1_img.to_filename(saveF)
    plotTitle = 'Corr coeffs, group, joint-ind contrast, right'
    savePlotF = folder + '/RsqMeanContrast1_img_real_right.svg'
    nPlotting.plot_stat_map(stat_map_img=saveF,
                            bg_img=backg_img,
                            threshold=0.10,
                            colorbar=True,
                            display_mode='ortho',
                            cut_coords=[55,-5,10],
                            annotate=True,
                            title=plotTitle,
                            output_file=savePlotF)

    # Joint stories, group
    RsqMeanContrast2 = part['RsqMeanContrast2']
    RsqMeanContrast2 = np.power(RsqMeanContrast2, 0.5)
    RsqMeanContrast2_img = masker.inverse_transform(RsqMeanContrast2)
    saveF = folder + '/RsqMeanContrast2_img_real.nii'
    RsqMeanContrast2_img.to_filename(saveF)
    plotTitle = 'Corr coeffs, group, ind-joint contrast, right'
    savePlotF = folder + '/RsqMeanContrast2_img_real_right.svg'
    nPlotting.plot_stat_map(stat_map_img=saveF,
                            bg_img=backg_img,
                            threshold=0.10,
                            colorbar=True,
                            display_mode='ortho',
                            cut_coords=[55,-5,10],
                            annotate=True,
                            title=plotTitle,
                            output_file=savePlotF)

return




def plotMeanRsqPerm():

    #%% File handling

    # get to right pwd
    folder = '/home/adamb/Desktop/testConfounds'
    os.chdir(folder)

    # specify background img, brainmask img, tr, verbosity
    backg_img = folder + '/mni_icbm152_t1_tal_nlin_asym_09c.nii'
    brainmask_img = folder + '/mni_asym09c_mask_resamp3x3.nii.gz'
    tr = 1.9
    masker_verbose = 2

    ## specify our target story
    #pairN = 2
    #runN = 1

    # exact mat files
    file = folder + '/permResults_First/permResults.mat'

    # load data
    part = sio.loadmat(file)


    #%% Prepare masker

    # create masker object using the 1x1x1 mm MNINlin mask
    masker = nInput.NiftiMasker(brainmask_img,
                                t_r=tr,
                                verbose=masker_verbose)
    # we need to fit the anatomical to it so we can use it with other data as well
    anatImg = masker.fit(backg_img)


    #%% Fit what we want and plot it

    # Joint stories, group
    RsqImg = part['RsqImg']
#    RsqMeanJoint = np.power(RsqMeanJoint, 0.5)
    RsqMeanImg = masker.inverse_transform(RsqImg)
    saveF = folder + '/permResults_First/RsqMeanImg_real.nii'
    RsqMeanImg.to_filename(saveF)
    plotTitle = 'R squared, group level, left'
    savePlotF = folder + '/permResults_First/RsqMeanImg_real_left.svg'
    nPlotting.plot_stat_map(stat_map_img=saveF,
                            bg_img=backg_img,
                            threshold=0.1,
                            colorbar=True,
                            display_mode='ortho',
                            cut_coords=[-55,-5,10],
                            annotate=True,
                            title=plotTitle,
                            output_file=savePlotF)

#%%
    return



def perm2nifti():
    #%% File handling

    # get to right pwd
    folder = '/home/adamb/Desktop/testConfounds'
    os.chdir(folder)
    
    # specify background img, brainmask img, tr, verbosity
    backg_img = folder + '/mni_icbm152_t1_tal_nlin_asym_09c.nii'
    brainmask_img = folder + '/mni_asym09c_mask_resamp3x3.nii.gz'
    fit_img = folder + '/results_First/singleEPIMNIsample.nii'
    tr = 1.9
    masker_verbose = 2

    ## specify our target story
    #pairN = 2
    #runN = 1

    # exact mat files
    file = folder + '/results_First/permResults.mat'

    # load data
    part = sio.loadmat(file)


    #%% Prepare masker

    # create masker object using the 3x3x3 mm MNINlin mask
    masker = nInput.NiftiMasker(brainmask_img,
                                t_r=tr,
                                verbose=masker_verbose)
    # we need to fit the anatomical to it so we can use it with other data as well
    sampleImg = masker.fit(fit_img)


    #%% Fit what we want and plot it

    # part 1
    RsqFDR = part['RsqImg']
    RsqFDRimg = masker.inverse_transform(RsqFDR)
    saveF = folder + '/results_First/permResults.nii.gz'
    RsqFDRimg.to_filename(saveF)

    return



