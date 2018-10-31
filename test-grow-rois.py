# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 12:47:08 2018

@author: onerva
"""

import functions
ROIInfoFile = '/media/onerva/KINGSTON/test-data/group_roi_mask-30-4mm_with_subcortl_and_cerebellum.mat'
savePath = '/media/onerva/KINGSTON/test-data/spherical-rois-test.nii'

ROICentroids, _, voxelCoordinates = functions.readROICentroids(ROIInfoFile, readVoxels=True)
ROIInfo = functions.growROIs(ROICentroids,voxelCoordinates)
functions.createNii(ROIInfo, savePath, imgSize=[45,54,45])
