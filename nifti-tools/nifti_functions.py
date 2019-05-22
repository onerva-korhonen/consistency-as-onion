# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 15:47:57 2018

@author: onerva

A set of functions for reading, writing and manipulating NIFTI files. These functions
largely utilize the nibabel library (http://nipy.org/nibabel/)
"""
import nibabel as nib
import numpy as np

def createNii(ROIInfo, savePath, imgSize=[45,54,45], affine=np.eye(4)):
    """
    Based on the given ROIInfo, creates a ROI mask (in .nii format) and saves 
    it to a given path. To construct the image of the NIFTI file, a 3D
    matrix of zeros is created and values of voxels belonging to a ROI are set to
    the index of this ROI.
    
    Parameters:
    -----------
    ROIInfo: dict, contains:
                   ROIMaps: list of ROISizes x 3 np.arrays, coordinates of voxels
                           belonging to each ROI. len(ROIMaps) = nROIs.
    savePath: str, path for saving the mask.
    imgSize: list of ints, dimensions of the image in the nii file. If saving a
             ROI mask created based on existing mask, set to the dimensions of the
             existing mask. (Default: [45,54,45]; the dimensions of the 4mm 
             Brainnetome mask)
    affine: np.array, an image coordination transformation (affine) matrix. (Default:
            an identity matrix)
             
    Returns:
    --------
    No direct output, saves the mask in NIFTI format to the given path
    """
    data = np.zeros(imgSize)
    ROIMaps = ROIInfo['ROIMaps']
    for i, ROI in enumerate(ROIMaps):
        if len(ROI.shape) == 1:
            data[ROI[0],ROI[1],ROI[2]] = i + 1
        else:
            for voxel in ROI:   
                data[voxel[0],voxel[1],voxel[2]] = i + 1
    img = nib.Nifti1Image(data,affine)     
    nib.save(img,savePath)


